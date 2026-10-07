# Utilisation du pipeline ARN

## Sources et validation

```bash
python scripts/run_pipeline.py --accession NR_000027.1 --output-dir outputs
python scripts/run_pipeline.py --sequence AUGGCUACG --label exemple --output-dir outputs_sequence
python scripts/run_pipeline.py --fasta data/examples/reference.fasta --mutant-fasta data/examples/mutant.fasta --output-dir outputs_fasta
```

Choisir exactement une source de référence : `--accession`, `--sequence`, `--fasta`
ou `--batch-fasta`. Pour une comparaison unique, ajouter `--mutant-sequence`
ou `--mutant-fasta`. Les substitutions sont positionnelles, sans alignement.

Les entrées sont converties en majuscules, les espaces ignorés, T transformé en U.
Seules A/C/G/T/U sont acceptées ; une base ambiguë, un gap ou une entrée vide produit
une erreur explicite avec sa position. Un FASTA unique nécessite un en-tête `>`
et une seule séquence. Son premier mot fournit l'identifiant ; `--label` le remplace.
La description et la provenance restent présentes dans le rapport.

Les noms de fichiers sont assainis pour Windows/Linux, sans modifier l'identifiant
biologique dans les données. Pour les chemins avec espaces, utiliser des guillemets.

## Cache NCBI

Le pipeline conserve FASTA normalisé et métadonnées dans `data/cache`.
`--cache-dir` choisit le répertoire ; `--refresh-cache` force un nouveau téléchargement.
Une entrée valide n'expire pas automatiquement. Une accession non versionnée peut
donc devenir obsolète ; la provenance conserve accession demandée/résolue, base,
URL, date UTC, SHA-256 et indicateur de lecture du cache.

Un cache corrompu est retéléchargé ; un rafraîchissement en échec conserve l'ancien
fichier mais fait échouer la commande. Les écritures sont atomiques. Un résumé NCBI
vide n'invalide pas une séquence valide ; la provenance reste disponible.
Les erreurs HTTP/connexion/délai sont propagées ; email et clé API ne sont pas
stockés dans le cache. Voir [les détails](releases/v0.4.1.md).

## Lots FASTA

```bash
python scripts/run_pipeline.py --batch-fasta data/examples/batch.fasta --output-dir outputs_batch
```

Le dossier de sortie doit être nouveau. Chaque entrée possède un dossier numéroté,
ses rapports et sa provenance. Les identifiants dupliqués n'écrasent aucune entrée.
Une erreur de séquence produit `error.json` et le traitement continue.
`batch_summary.json` récapitule succès, erreurs et chemins relatifs. Code de sortie
0 en cas de succès complet, 2 si une entrée échoue. Les erreurs de cadrage FASTA
ou d'accès disque interrompent le lot. Ce mode ne compare pas de mutants.

Le calcul est séquentiel et conserve les entrées en mémoire ; le prédicteur
Nussinov a un coût cubique en temps et quadratique en mémoire.

## Résultats et réglages

```text
outputs/
  identifiant.json
  identifiant.md
  identifiant_hotspots.csv
  identifiant/visualization/
    visualization_manifest.json
    secondary_structures.json
    wt_structure.pdb
    mutant_structure.pdb
    wt_scores.json
    mutant_scores.json
    impact_summary.json
    hotspots.csv
```

Les fichiers mutants n'existent que si une comparaison est demandée.
`--top-k` (10), `--min-abs-delta` / `--delta-threshold` (0,1) et
`--base-pair-threshold` (0,2) règlent le résumé des variations. Les hotspots commencent
à 1 ; certains champs JSON historiques commencent à 0 : voir [les coordonnées](coordinates.md).

## Visualisation

```bash
python scripts/run_demo.py
python scripts/view_secondary.py --manifest outputs_demo/offline_demo/visualization/visualization_manifest.json
python scripts/view_secondary.py --manifest outputs_demo/offline_demo/visualization/visualization_manifest.json --export outputs_demo/secondary.png --export-only
python scripts/view_3d.py --manifest outputs_demo/offline_demo/visualization/visualization_manifest.json
```

La vue secondaire représente les paires prédites, leurs pertes/gains et les deltas.
Les PDB du pipeline ARN suivent une géométrie illustrative dépendant uniquement de
la longueur, sans prédire le repliement. Ils sont distincts de la reconstruction
Hi-C par contacts introduite en v0.4.6. Voir [le guide des viewers](visualization_3d.md).
