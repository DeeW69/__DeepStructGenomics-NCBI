"""Exercise viewer argument handling without importing Tkinter or VTK viewers."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
from types import ModuleType
from unittest.mock import Mock

import pytest


@pytest.fixture
def viewer_cli(monkeypatch):
    constructors = {}
    exports = {
        "tk_vtk_minimal": ("MinimalVTKViewer",),
        "tk_vtk_molecule": ("MoleculeViewer",),
        "tk_vtk_overlay": ("BasePairOptions", "OverlayViewer", "OverlayDeltaViewer"),
    }
    for module_name, class_names in exports.items():
        full_name = f"deepstructgenomics.visualization.{module_name}"
        module = ModuleType(full_name)
        for class_name in class_names:
            constructor = Mock(name=class_name)
            setattr(module, class_name, constructor)
            constructors[class_name] = constructor
        monkeypatch.setitem(sys.modules, full_name, module)

    script_path = Path(__file__).resolve().parents[1] / "scripts" / "view_3d.py"
    spec = importlib.util.spec_from_file_location("viewer_cli_under_test", script_path)
    cli = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cli)
    return cli, constructors


@pytest.mark.parametrize(
    "mode,arguments,class_name,expected_metric",
    [
        ("minimal", [], "MinimalVTKViewer", None),
        ("molecule", ["--structure", "molecule.pdb"], "MoleculeViewer", None),
        (
            "overlay",
            ["--wt", "wt.pdb", "--mut", "mut.pdb", "--scores", "scores.json"],
            "OverlayViewer",
            "score",
        ),
        (
            "overlay-delta",
            ["--wt", "wt.pdb", "--mut", "mut.pdb", "--wt-scores", "wt.json", "--mut-scores", "mut.json"],
            "OverlayDeltaViewer",
            "delta",
        ),
    ],
)
def test_default_link_metric_allows_each_viewer_mode(
    viewer_cli, monkeypatch, mode, arguments, class_name, expected_metric
):
    cli, constructors = viewer_cli
    monkeypatch.setattr(sys, "argv", ["view_3d.py", "--mode", mode, *arguments])

    cli.main()

    constructor = constructors[class_name]
    constructor.assert_called_once()
    constructor.return_value.start.assert_called_once()
    if mode == "minimal":
        constructor.assert_called_once_with()
    elif mode == "molecule":
        constructor.assert_called_once_with("molecule.pdb", show_backbone=False)
    else:
        assert constructor.call_args.kwargs["link_metric"] == expected_metric


@pytest.mark.parametrize("mode", ["minimal", "molecule", "overlay"])
def test_explicit_delta_metric_still_requires_delta_mode(viewer_cli, monkeypatch, mode):
    cli, constructors = viewer_cli
    monkeypatch.setattr(sys, "argv", ["view_3d.py", "--mode", mode, "--link-metric", "delta"])

    with pytest.raises(SystemExit, match="--link-metric delta requiert --mode overlay-delta"):
        cli.main()

    for constructor in constructors.values():
        constructor.assert_not_called()
