![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# DeepStructGenomics-NCBI

DeepStructGenomics-NCBI est un cadre de recherche et d'ingénierie logicielle dédié à l'étude structurale de l'ARN et de la régulation génomique à partir des données publiques du NCBI. L'objectif est de relier systématiquement séquences biologiques, structures (ARN et architecture 3D) et impact fonctionnel potentiel, afin de proposer des analyses mécanistiques exploitables en recherche biomédicale.

## Objectifs scientifiques

- **Pipeline séquence → structure → impact** : ingestion de séquences ou d'identifiants NCBI, annotation, modélisation structurale (surtout ARN) puis génération d'indicateurs fonctionnels.
- **Interprétation des variants** : comparaison de séquences référence/mutant pour estimer les effets structuraux ou régulateurs.
- **Génomique 3D** (à venir) : synthèse de données de contacts (Hi-C et apparentées) pour relier organisation spatiale du génome et régulation.

## Architecture de dépôt

```
DeepStructGenomics-NCBI/
|-- README.md
|-- requirements.txt
|-- pyproject.toml
|-- docs/
|   |-- visualization_3d.md
|-- data/
|   |-- examples/
|-- deepstructgenomics/
|   |-- __init__.py
|   |-- config.py
|   |-- pipeline.py
|   |-- data_sources/
|   |   |-- __init__.py
|   |   |-- ncbi_client.py
|   |-- rna/
|   |   |-- __init__.py
|   |   |-- secondary_structure.py
|   |-- variants/
|   |   |-- __init__.py
|   |   |-- variant_analysis.py
|   |-- reporting/
|   |   |-- __init__.py
|   |   |-- report_generator.py
|   |-- visualization/
|       |-- __init__.py
|       |-- io_structures.py
|       |-- score_mapping.py
|       |-- tk_vtk_minimal.py
|       |-- tk_vtk_molecule.py
|       |-- tk_vtk_overlay.py
|-- scripts/
|   |-- run_pipeline.py
|   |-- view_3d.py
```

- `data_sources` couvre l'interfaçage NCBI (séquences, annotations, variants).
- `rna` contient les modèles de structures secondaires et l'extraction de motifs.
- `variants` fournit les comparaisons référence/mutant (impact structural).
- `reporting` regroupe la génération de sorties JSON et Markdown.
- `visualization` implémente la lecture des structures, la conversion de scores et les viewers Tkinter + VTK.
- `scripts/run_pipeline.py` expose la CLI principale, `scripts/view_3d.py` lance les viewers 3D.

## Pipeline minimal (V1)

Entrées possibles :

1. Identifiant NCBI (nucléotidique, par ex. `NR_000027.1`).
2. Séquence ARN fournie directement, accompagnée d'un nom ou d'un identifiant libre.

Sorties :

- Fichier JSON structuré incluant métadonnées, annotations simples, score de stabilité, paires de bases prédites.
- Rapport Markdown synthétisant les mêmes informations et intégrant une représentation ASCII de la structure secondaire.
- (optionnel) Dossier `visualization/` avec structures pseudo-3D (WT/mutant), scores normalisés et manifeste.

### Étapes principales

1. **Récupération / normalisation** : téléchargement via E-utilities (NCBI) ou validation de la séquence fournie.
2. **Annotation simple** : calcul longueur, composition, GC%, heuristiques de fonctionnalités.
3. **Prédiction de structure secondaire** : implémentation interne d'un algorithme type Nussinov avec matrice énergie simplifiée.
4. **Analyse variant (optionnelle V1)** : si une séquence mutante est fournie, comparaison des structures prédites et estimation heuristique de l'impact.
5. **Reporting et visualisation** : sérialisation JSON/Markdown + génération de structures pseudo-3D et scores dérivés.

## Format de sortie

```text
outputs/
  NR_000027.1.json
  NR_000027.1.md
  NR_000027.1/
    visualization/
      wt_structure.pdb
      mutant_structure.pdb       # si variant fourni
      mutant_scores.json         # scores [0..1] par résidu / position
      visualization_manifest.json
```

Les structures 3D générées sont des projections contraintes (une sphère par nucléotide) destinées à comparer des états WT/mutant. Elles ne représentent pas une conformation atomique unique. Les scores sont dérivés des appariements (structure secondaire) et servent à visualiser des gradients relatifs.

## Installation rapide

```bash
python -m venv .venv
source .venv/bin/activate  # sous Windows: .\.venv\Scripts\activate
pip install -r requirements.txt
```

## Utilisation du pipeline

```bash
python scripts/run_pipeline.py --accession NR_000027.1 --output-dir outputs
# ou
python scripts/run_pipeline.py --sequence AUGGCUACG --label test_seq --output-dir outputs
```

Les rapports `.json` et `.md` sont écrits dans le dossier de sortie, tandis que les artefacts de visualisation sont placés dans `outputs/<identifiant>/visualization/`.

## Visualisation 3D contrainte

Le viewer desktop est accessible via :

```bash
python scripts/view_3d.py --mode minimal
python scripts/view_3d.py --mode molecule --structure path/to/model.pdb
python scripts/view_3d.py --mode overlay --wt path/to/wt.pdb --mut path/to/mut.pdb --scores path/to/mutant_scores.json
```

Consultez `docs/visualization_3d.md` pour les détails d'installation, les formats attendus (`mutant_scores.json`, `visualization_manifest.json`) et les limites scientifiques (structure secondaire = socle, projection 3D = contrainte, détails tertiaires = exploratoires).

## Principes directeurs

- **Reproductibilité** : configuration explicite, dépendances verrouillées, logs détaillés.
- **Traçabilité** : chaque étape du pipeline enregistre ses paramètres et résultats intermédiaires (incluant les artefacts de visualisation).
- **Modularité** : les modules peuvent évoluer indépendamment pour intégrer de nouvelles méthodes (predictors ARN, intégration Hi-C).
- **Explicabilité** : prioriser des modèles interprétables, exposer les hypothèses biologiques et les incertitudes des prédictions.
- **Open science** : privilégier données publiques, formats ouverts et documentation accessible.
