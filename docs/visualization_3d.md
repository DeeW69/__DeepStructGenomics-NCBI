# Visualisation 3D (Tkinter + VTK)

Ce module fournit une visualisation de laboratoire pour inspecter les projections 3D contraintes utilisees par DeepStructGenomics. Il s'agit d'une representation fonctionnelle derivee de la structure secondaire ARN validee : l'objectif est de raisonner sur les voisinages, l'accessibilite et les variations relatives WT/mutant, pas de decrire une conformation atomiquement exacte.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .\.venv\Scripts\activate
pip install -r requirements.txt
# ou au minimum :
pip install vtk biopython numpy
```

Ces dependances sont necessaires pour charger des structures PDB/mmCIF (BioPython), manipuler des scalaires (NumPy) et afficher les glyphes 3D (VTK).

## CLI `scripts/view_3d.py`

```bash
# Viewer minimal : camera + acteur de demonstration
python scripts/view_3d.py --mode minimal

# Viewer molecule : chargement direct d'un PDB/mmCIF
python scripts/view_3d.py --mode molecule --structure data/examples/mini_helix.pdb

# Overlay WT/mutant : necessite deux structures et un fichier de scores JSON
python scripts/view_3d.py --mode overlay \
  --wt outputs/<run_id>/visualization/wt_structure.pdb \
  --mut outputs/<run_id>/visualization/mutant_structure.pdb \
  --scores outputs/<run_id>/visualization/mutant_scores.json
```

Chaque fenetre est une application Tkinter autonome. Les interactions standards VTK sont disponibles (rotation, zoom molette, translation + clic droit).

## Flux pipeline

Lors d'un `python scripts/run_pipeline.py ...`, le pipeline ecrit desormais :

```text
outputs/
  NR_000027.1.json
  NR_000027.1.md
  NR_000027.1/
    visualization/
      wt_structure.pdb
      mutant_structure.pdb        # si sequence mutante fournie
      mutant_scores.json          # scores [0..1] par residu / position
      visualization_manifest.json
```

- `wt_structure.pdb` et `mutant_structure.pdb` sont des projections helicoidales contraintes (1 atome/glyphe par nucleotide).
- `mutant_scores.json` encode les scores de stabilite/accessibilite derives de la structure secondaire. Les cles peuvent etre soit `position_scores` (index 1..N), soit `residue_scores` (cle `chaine:resseq:icode`). Tout atome absent recoit un score 0.
- `visualization_manifest.json` fournit la tracabilite : identifiant NCBI, source, date, parametres du modele, fichiers produits.

### Exemple de `mutant_scores.json`

```json
{
  "version": 1,
  "metadata": {
    "context": "mutant",
    "reference_identifier": "NR_000027.1"
  },
  "residue_scores": {
    "A:1:": 0.65,
    "A:2:": 0.91
  },
  "position_scores": {
    "1": 0.65,
    "2": 0.91
  }
}
```

## Limites scientifiques & hierarchie de confiance

- **Structure secondaire = socle (niveau eleve)** : les scores sont calcules a partir des appariements dot-bracket valides. Toute interpretation doit partir de ce niveau.
- **Projection 3D contrainte (niveau intermediaire)** : la geometrie est generee par une helice parametrique pour visualiser des voisinages et contraintes spatiales. Ce n'est **pas** une structure atomique physique.
- **Details tertiaires (niveau exploratoire)** : couleurs / glyphes servent a comparer WT vs mutant. Les conclusions doivent rester comparatives (ou les scores evoluent) plutot qu'absolues.

Les metadonnees fournies (identifiant NCBI, parametres de generation, datation UTC) permettent de retracer les hypotheses et de reproduire la visualisation. Toute utilisation biomedicale doit tenir compte du contexte cellulaire, des partenaires moleculaires et des incertitudes sur l'environnement.
