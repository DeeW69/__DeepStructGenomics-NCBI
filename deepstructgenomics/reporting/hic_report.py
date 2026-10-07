"""Portable tables and static research figures for advanced Hi-C analysis."""
import csv
import json
from pathlib import Path

import numpy as np


def export_hic_analysis(report, output_dir, *, plot=True):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "hic_analysis.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    region = report["region"]
    start, size, chromosome = region["start"], region["bin_size"], region["chromosome"]

    def table(name, columns, rows):
        with (output_dir / name).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    matrix = np.asarray(report["balanced_matrix"])
    table("balanced_contacts.tsv", ["chrom1", "start1", "chrom2", "start2", "weight"],
          ({"chrom1": chromosome, "start1": start + i * size, "chrom2": chromosome,
            "start2": start + j * size, "weight": float(matrix[i, j])}
           for i in range(len(matrix)) for j in range(i + 1, len(matrix))))
    if report["status"] != "ok":
        (output_dir / "hic_analysis.md").write_text(
            "# Hi-C : équilibrage non convergé\n\nAucun test statistique ni reconstruction exécuté.\n",
            encoding="utf-8")
        return
    loops, domains, reconstruction = report["loops"], report["domains"], report["reconstruction"]
    table("loop_tests.tsv", ["chromosome", "start1", "start2", "count", "expected", "dispersion", "enrichment", "p_value", "q_value", "significant"],
          ({"chromosome": chromosome, "start1": start + row["bin1"] * size, "start2": start + row["bin2"] * size,
            **{key: row[key] for key in ("count", "expected", "dispersion", "enrichment", "p_value", "q_value", "significant")}}
           for row in loops["tests"]))
    table("boundary_tests.tsv", ["chromosome", "position", "insulation", "p_value", "q_value", "selected"],
          ({"chromosome": chromosome, "position": start + row["boundary_bin"] * size,
            **{key: row[key] for key in ("insulation", "p_value", "q_value", "selected")}}
           for row in domains["tests"]))
    table("coordinates.tsv", ["chromosome", "start", "end", "x", "y", "z"],
          ({"chromosome": chromosome, "start": start + row["bin"] * size, "end": start + (row["bin"] + 1) * size,
            **{key: row[key] for key in ("x", "y", "z")}} for row in reconstruction["coordinates"]))
    with (output_dir / "candidate_domains.bed").open("w", encoding="utf-8") as stream:
        for number, interval in enumerate(domains["candidate_domains"], 1):
            stream.write(f"{chromosome}\t{start + interval['start_bin'] * size}\t{start + interval['end_bin'] * size}\tcandidate_{number}\n")
    significant = sum(row["significant"] for row in loops["tests"])
    lines = ["# Analyse Hi-C exploratoire", "",
             f"Région : {chromosome}:{start}-{region['end']} ({report['assembly']}), bins de {size} bp.",
             f"Équilibrage convergé en {report['normalization']['iterations']} itérations ; "
             f"erreur relative maximale {report['normalization']['max_relative_row_error']:.3g}.",
             f"Boucles candidates : {significant} / {loops['tested_pixels']} pixels testés.",
             f"Frontières retenues : {len(domains['selected_boundaries'])}.",
             f"Reconstruction : {reconstruction['status']} ; stress : {reconstruction.get('stress', 'indisponible')}.",
             "", "Les q-values dépendent des modèles nuls décrits dans le JSON. Les intervalles",
             "sont des domaines candidats ; les coordonnées sont inférées en unités arbitraires.",
             "Aucune validation biologique ni calibration physique n'est impliquée.", ""]
    (output_dir / "hic_analysis.md").write_text("\n".join(lines), encoding="utf-8")
    if plot:
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        figure = Figure(figsize=(11, 4.5), constrained_layout=True)
        FigureCanvasAgg(figure)
        axis = figure.add_subplot(121)
        displayed = matrix.copy()
        masked = report["normalization"]["masked_bins"]
        displayed[masked, :] = np.nan
        displayed[:, masked] = np.nan
        artist = axis.imshow(displayed, origin="lower", cmap="magma")
        figure.colorbar(artist, ax=axis, label="Balanced contact weight")
        for boundary in domains["selected_boundaries"]:
            axis.axvline(boundary - 0.5, color="cyan", linewidth=0.8)
            axis.axhline(boundary - 0.5, color="cyan", linewidth=0.8)
        axis.set(title="Balanced cis contacts", xlabel="Bin index (0-based)", ylabel="Bin index (0-based)")
        axis3d = figure.add_subplot(122, projection="3d")
        coords = reconstruction["coordinates"]
        if coords:
            xyz = np.array([[row[key] for key in ("x", "y", "z")] for row in coords])
            axis3d.scatter(*xyz.T, c=[row["bin"] for row in coords], cmap="viridis", s=18)
            for i in range(len(coords) - 1):
                if coords[i + 1]["bin"] == coords[i]["bin"] + 1:
                    axis3d.plot(*xyz[i:i + 2].T, color="gray", linewidth=0.6)
        axis3d.set(title=f"Inferred MDS — {reconstruction['status']}", xlabel="x (a.u.)", ylabel="y (a.u.)", zlabel="z (a.u.)")
        figure.savefig(output_dir / "hic_analysis.png", dpi=150)
