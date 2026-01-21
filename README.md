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
├── README.md
├── requirements.txt
├── pyproject.toml
├── data/
│   └── examples/
├── deepstructgenomics/
│   ├── __init__.py
│   ├── config.py
│   ├── pipeline.py
│   ├── data_sources/
│   │   ├── __init__.py
│   │   └── ncbi_client.py
│   ├── rna/
│   │   ├── __init__.py
│   │   └── secondary_structure.py
│   ├── variants/
│   │   ├── __init__.py
│   │   └── variant_analysis.py
│   └── reporting/
│       ├── __init__.py
│       └── report_generator.py
└── scripts/
    └── run_pipeline.py
```

- `data_sources` couvre l'interfaçage NCBI (séquences, annotations, variants).
- `rna` contient les modèles de structures secondaires et l'extraction de motifs.
- `variants` fournit les comparaisons simple référence/mutant (impact structural).
- `reporting` regroupe la génération de sorties JSON et Markdown.
- `pipeline.py` orchestre l'enchaînement reproductible des modules.
- `scripts/run_pipeline.py` expose une CLI minimale pour la V1.

## Pipeline minimal (V1)

Entrées possibles :

1. Identifiant NCBI (nucléotidique, par ex. `NR_000027.1`).
2. Séquence ARN fournie directement, accompagnée d'un nom ou d'un identifiant libre.

Sorties :

- Fichier JSON structuré incluant métadonnées, annotations simples, score de stabilité, paires de bases prédites.
- Rapport Markdown synthétisant les mêmes informations et intégrant une représentation ASCII de la structure secondaire.

### Étapes principales

1. **Récupération / normalisation** : téléchargement via E-utilities (NCBI) ou validation de la séquence fournie.
2. **Annotation simple** : calcul longueur, composition, GC%, heuristiques de fonctionnalités.
3. **Prédiction de structure secondaire** : implémentation interne d'un algorithme type Nussinov avec matrice énergie simplifiée pour produire des paires base-base et des motifs élémentaires.
4. **Analyse variant (optionnelle V1)** : si une séquence mutante est fournie, comparaison des structures prédites et estimation heuristique de l'impact.
5. **Reporting** : sérialisation JSON + gabarit Markdown.

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

Les résultats sont écrits sous forme `*.json` et `*.md` dans le dossier de sortie. Les fichiers sont autoportants pour l'intégration en pipelines aval (CI/CD, notebooks d'analyse, bases de connaissances).

## Principes directeurs

- **Reproductibilité** : configuration explicite, dépendances verrouillées, logs détaillés.
- **Traçabilité** : chaque étape du pipeline enregistre ses paramètres et résultats intermédiaires.
- **Modularité** : les modules peuvent évoluer indépendamment pour intégrer de nouvelles méthodes (ex. prédicteurs ARN plus précis, intégration Hi-C).
- **Explicabilité** : prioriser des modèles interprétables, exposer les hypothèses biologiques et les incertitudes des prédictions.
- **Open science** : privilégier données publiques, formats ouverts et documentation accessible.
