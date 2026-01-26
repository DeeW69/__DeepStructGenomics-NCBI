# Visualisation 3D (VTK standalone)

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

# Overlay delta (comparatif WT -> mutant)
python scripts/view_3d.py --mode overlay-delta \
  --wt outputs/<run_id>/visualization/wt_structure.pdb \
  --mut outputs/<run_id>/visualization/mutant_structure.pdb \
  --wt-scores outputs/<run_id>/visualization/wt_scores.json \
  --mut-scores outputs/<run_id>/visualization/mutant_scores.json
```

Chaque fenetre est une application VTK standalone (interactor natif). Les interactions standards VTK sont disponibles (rotation, zoom molette, translation + clic droit).

Options utiles :

- `--no-backbone` desactive le backbone discret pour overlay / overlay-delta (active par defaut).
- `--backbone` force le backbone en mode molecule.
- `--export chemin.png` capture un screenshot des l'ouverture.
- `--export-only` sauvegarde le screenshot puis quitte sans lancer l'interacteur (pratique en CI).

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
      wt_scores.json              # scores [0..1] pour la reference
      mutant_scores.json          # scores [0..1] pour le mutant
      visualization_manifest.json
```

- `wt_structure.pdb` et `mutant_structure.pdb` sont des projections helicoidales contraintes (1 atome/glyphe par nucleotide).
- `wt_scores.json` et `mutant_scores.json` encodent les scores derives de la structure secondaire (cles `position_scores` et/ou `residue_scores`). Ils sont utilises pour le mode overlay (`mutant_scores.json`) et pour le delta (les deux fichiers).
- `visualization_manifest.json` fournit la tracabilite : identifiant NCBI, source, date, parametres du modele, fichiers produits, statistiques (`score_statistics_*`) ainsi que le mode d'etirement applique (note : la LUT etend le contraste visuel mais ne modifie pas les donnees).

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
- **Details tertiaires (niveau exploratoire)** : couleurs / glyphes servent a comparer WT vs mutant. Les conclusions doivent rester comparatives (delta = variation relative) plutot qu'absolues.

Les metadonnees fournies (identifiant NCBI, parametres de generation, datation UTC) permettent de retracer les hypotheses et de reproduire la visualisation. Toute utilisation biomedicale doit tenir compte du contexte cellulaire, des partenaires moleculaires et des incertitudes sur l'environnement.

## Debug / diagnostics

Mettre `DSG_VIZ_DEBUG=1` dans l'environnement avant `python scripts/view_3d.py ...` permet d'afficher des statistiques (couverture des scores, range utilise pour la LUT) lors du rendu overlay. Utile pour verifier l'alignement des indices residus ↔ atomes, sans polluer la sortie en mode normal.

## Mode overlay-delta

Le mode `overlay-delta` colore les points du mutant selon `delta = score_mutant - score_wt` :

- bleu = diminution / destabilisation relative,
- rouge = augmentation / stabilisation relative,
- gris clair autour de 0 (impact neutre).

La LUT est centree sur 0 et etiree de maniere symetrique (`range_mode: "dynamic_symmetric"` dans le manifeste). Les scores bruts ne sont pas modifies ; seul l'affichage ajuste la plage pour rendre le gradient lisible. Utiliser ce mode uniquement lorsque les deux fichiers `wt_scores.json` et `mutant_scores.json` proviennent du meme run afin de garantir la coherence des index.

## Legende & tooltip

Les modes overlay affichent une legende textuelle rappelant que WT est rendu en gris translucide tandis que le mutant est colore selon son score (0..1) ou `Δ = mutant - WT`. La barre graduée indique la plage reelement etiree (symetrique autour de 0 pour le delta).  
Pour activer un tooltip detaille (position, residu, valeurs WT/MUT/Δ et range de clamp), definir `DSG_VIZ_DEBUG=1` avant le lancement : un clic/survol sur un point coloré afﬁchera ces informations. Sans cette variable, aucun tooltip n'est installe afin de garder la fenetre silencieuse.

## Backbone discret

Overlay et overlay-delta ajoutent un backbone/tube gris translucide pour suivre la continuite des pseudo-atomes. Il peut etre desactive via `--no-backbone`. Le mode molecule peut afficher le meme tube avec `--backbone`.

## Export PNG

Toutes les commandes acceptent `--export chemin.png` pour capturer une image. Avec `--export-only`, la capture est effectuee puis l'application se ferme sans lancer l'interacteur.

```powershell
python scripts/view_3d.py --mode overlay `
  --wt outputs/runA/visualization/wt_structure.pdb `
  --mut outputs/runA/visualization/mutant_structure.pdb `
  --scores outputs/runA/visualization/mutant_scores.json `
  --export outputs/runA/visualization/overlay.png `
  --export-only
```
