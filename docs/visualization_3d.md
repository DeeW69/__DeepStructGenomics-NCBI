# Visualisation 3D (VTK standalone)

Ce module affiche les PDB/mmCIF et les scores du pipeline. Les PDB générés
utilisent une hélice paramétrique dépendant seulement de la longueur de la
séquence, sans repliement moléculaire prédit. Des séquences de même longueur
ont les mêmes coordonnées ; les distances affichées ne mesurent donc pas un
déplacement biologique. Les scores et les connecteurs proviennent des
appariements secondaires prédits par l'heuristique.

Pour une comparaison lisible des appariements, privilégier
`scripts/view_secondary.py --manifest chemin/visualization_manifest.json`
après régénération des résultats sur `main` ; voir le [guide de démo](demo.md).
Ce viewer ouvre désormais les tiges/boucles 2D par défaut ; `--view arcs` conserve
le diagramme linéaire et les deltas. Dans le desktop, l'onglet **3D schématique**
reste complémentaire à **Structure secondaire**.

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
# Chargement automatique depuis un manifeste (recommandé) :
python scripts/view_3d.py --manifest outputs/<run_id>/visualization/visualization_manifest.json

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

- `--manifest fichier.json` résout automatiquement les structures, scores et connecteurs, et détecte le mode (`overlay-delta`, `overlay` ou `molecule`).
- `--no-backbone` desactive le backbone discret pour overlay / overlay-delta (active par defaut).
- `--backbone` force le backbone en mode molecule.
- `--export chemin.png` capture un screenshot des l'ouverture.
- `--export-only` sauvegarde le screenshot puis quitte sans lancer l'interacteur (pratique en CI).
- `--export-view preset` force un preset camera pour le PNG (`auto`, `iso`, `top`, `side`, `front`).
- `--export-hide-ui` cache temporairement les legendes/barres/statuts pendant l'export (affichées à nouveau en interactif).
- `--export-background {dark,white,transparent}` choisit le fond du PNG (transparent active les canaux alpha).
- `--tooltips` installe un tooltip compact meme si `DSG_VIZ_DEBUG` n'est pas defini.
- `--tooltips-verbose` force la version detaillee du tooltip (sans logs) meme sans `DSG_VIZ_DEBUG`.
- `--export-scale N` multiplie la resolution du PNG exporte (`N=2` double largeur/hauteur).
- `--links` active l'affichage de segments WT->MUT (correspondance positionnelle).
- `--no-links` force la desactivation des liens meme si `--links` est present dans un script appelant.
- `--link-threshold T` (ex: `0.2`) n'affiche que les liens dont `|score|` (overlay) ou `|delta|` (overlay-delta) depasse `T`.
- `--link-metric {score,delta,distance}` change la metrique utilisee pour colorer les liens (`score` par defaut en mode overlay, `delta` en mode overlay-delta ; l'option `delta` n'est disponible qu'avec `--mode overlay-delta`).
- `--link-min-distance` / `--link-max-distance` filtrent les liens selon la distance WT→MUT en unit��s VTK (utile avec `--link-metric distance`).
- `--base-pairs-wt` / `--base-pairs-mut` pointent vers des JSON `{"base_pairs": [[i,j], ...]}` (indices 1-based) de la structure secondaire WT et/ou MUT.
- `--base-pairs-mode {lines,arcs}` change l'aspect des connecteurs (segments droits ou arcs discrets).
- `--base-pairs-threshold <float>` applique (mode overlay-delta) un filtre sur `max(|mut-wt|)` des positions i/j ; utile pour ne visualiser que les appariements les plus impactés.
- `--no-base-pairs` force la désactivation des connecteurs, même si des JSON sont fournis par défaut dans un script.

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
      impact_summary.json         # nouveau resume WT->MUT (toujours present)
      visualization_manifest.json
```

- `wt_structure.pdb` et `mutant_structure.pdb` sont des projections hélicoïdales illustratives (1 atome/glyphe par nucléotide), déterminées par la longueur.
- `wt_scores.json` et `mutant_scores.json` encodent les scores derives de la structure secondaire (cles `position_scores` et/ou `residue_scores`). Ils sont utilises pour le mode overlay (`mutant_scores.json`) et pour le delta (les deux fichiers).
- `visualization_manifest.json` fournit la tracabilite : identifiant NCBI, source, date, parametres du modele, fichiers produits, statistiques (`score_statistics_*`, `displacement_statistics`, etc.) ainsi que le mode d'etirement applique (note : la LUT etend le contraste visuel mais ne modifie pas les donnees). Le manifest reference maintenant explicitement `impact_summary_file`, ce qui permet aux outils en aval de retrouver facilement le resume d'impact.
- `impact_summary.json` contient un resume WT->MUT nullement lie a VTK : statistiques globales (min/max/moyenne/medianes des deltas, ratios positifs/negatifs), "hotspots" (top positions triees par `|delta|` au-dessus d'un seuil configurable) et, si les paires de bases sont connues, un comptage des connecteurs touches (`max(|delta(i)|, |delta(j)|) >= seuil`). Lorsque la sequence mutante n'est pas fournie, ce fichier existe quand meme avec des zeros et une note explicite (`"mutant not provided"`). Ce fichier est concu pour etre digeste rapidement dans des notebooks / scripts d'analyse.

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

Consulter le [guide de numérotation des positions](coordinates.md) pour la table
de correspondance complète entre numéros de résidus PDB, fichiers de scores,
connecteurs de paires et rapports.

## Limites scientifiques & hierarchie de confiance

- **Structure secondaire = socle (niveau eleve)** : les scores sont calcules a partir des appariements dot-bracket valides. Toute interpretation doit partir de ce niveau.
- **Projection 3D illustrative** : l'hélice paramétrique ne dépend pas des appariements et ne permet pas de déduire voisinage moléculaire, accessibilité ou contraintes spatiales physiques.
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

## Legende, status text et tooltips

Les modes overlay et overlay-delta affichent maintenant un scalar bar (meme LUT que les glyphes) et un status text compact dans le coin inferieur gauche. Ce texte liste le mode actif, les fichiers de scores utilises, la plage scalaire appliquee ainsi que la plage de clamp, rappelle la legende (WT gris translucide, mutant colore) et termine par les raccourcis `Keys: R reset | W WT | M MUT | B backbone | L links | P base pairs`. En mode overlay-delta, une ligne supplementaire precise `Delta: 0=neutral | blue=decrease red=increase`. Quand `--links` est actif une ligne ASCII ajoute `Links: WT -> MUT (same sequence position), color follows active LUT`. Si les paires de bases sont activees, une ligne `BasePairs: on/off | WT=<n> MUT=<n> | mode=<lines/arcs> | threshold=<value> | toggle=P` resume l'etat du rendu et rappelle le raccourci `P`.

Par defaut les tooltips restent desactives, sauf si `DSG_VIZ_DEBUG=1`. Utiliser `--tooltips` pour activer la version compacte (position + residu + valeur affichee). `--tooltips-verbose` ou `DSG_VIZ_DEBUG=1` affichent la version detaillee (valeurs WT/MUT/delta + ranges). Les traces console additionnelles ne sont affichees que lorsque `DSG_VIZ_DEBUG=1`.

## Liens WT→MUT (optionnels)

`--links` affiche des segments reliant chaque point WT a son homologue mutant (meme position de sequence ou index 1-based en repli). Les couleurs suivent par defaut la LUT active (score mutant pour overlay, delta pour overlay-delta) mais on peut choisir `--link-metric distance` pour colorer en fonction du deplacement WT→MUT (`--link-metric delta` est reserve au mode overlay-delta). Dans le mode distance, la barre de statut ajoute la ligne `links metric: distance (min=..., max=...)` et les tooltips (si actifs) affichent la distance locale. Ces segments sont purement informatifs: ils representent une correspondance positionnelle, pas une liaison chimique.

La touche `L` permet de masquer/afficher les liens a chaud. Plusieurs filtres sont disponibles:

- `--link-threshold` conserve uniquement les liens dont `|score|` ou `|delta|` depasse le seuil.
- `--link-min-distance` / `--link-max-distance` filtrent par distance WT→MUT (utile pour isoler les deplacements significatifs, meme si la metrique de couleur reste le score/delta).

Exemple pour visualiser les liens par distance:

```powershell
python scripts/view_3d.py --mode overlay `
  --wt outputs/.../wt_structure.pdb `
  --mut outputs/.../mutant_structure.pdb `
  --scores outputs/.../mutant_scores.json `
  --links `
  --link-metric distance `
  --link-min-distance 2.0
```

## Paires de bases (secondary structure)

Les connecteurs de paires de bases visualisent les appariements secondaires
(indices depuis 1), pas des liaisons chimiques 3D. Le pipeline exporte désormais
`wt_base_pairs.json` et `mut_base_pairs.json`, détectés automatiquement avec
`--manifest`. Avec les anciens résultats ou sans manifeste, fournir ces fichiers
via `--base-pairs-wt` et `--base-pairs-mut`. `--no-base-pairs` masque les connecteurs.
Le format attendu est minimal :

```json
{
  "base_pairs": [[5, 28], [6, 27], [10, 15]]
}
```

Les indices sont automatiquement tries (`[i,j]` devient `(min,max)`), deduples puis tries pour garantir une liste stable. Les connecteurs sont traces en gris clair (`--base-pairs-mode lines` ou `arcs`), et peuvent etre masques/aﬃches a chaud via la touche `P`. Quand `--base-pairs-threshold > 0` est utilise en mode `overlay-delta`, le moteur filtre les paires selon `max(|mut-wt|)` calcule sur les positions i/j (scores issus des fichiers `--wt-scores`/`--mut-scores`). Dans ce cas, si toutes les valeurs sont disponibles, les segments heritent de la LUT delta pour souligner les appariements fortement impactes ; sinon le filtre retombe automatiquement en mode conservateur (pas de suppression).

Exemples :

```powershell
python scripts/view_3d.py --mode overlay `
  --wt outputs/NR_000027.1/visualization/wt_structure.pdb `
  --mut outputs/NR_000027.1/visualization/mutant_structure.pdb `
  --scores outputs/NR_000027.1/visualization/mutant_scores.json `
  --base-pairs-wt outputs/NR_000027.1/visualization/wt_base_pairs.json `
  --base-pairs-mut outputs/NR_000027.1/visualization/mut_base_pairs.json

python scripts/view_3d.py --mode overlay-delta `
  --wt outputs/NR_000027.1/visualization/wt_structure.pdb `
  --mut outputs/NR_000027.1/visualization/mutant_structure.pdb `
  --wt-scores outputs/NR_000027.1/visualization/wt_scores.json `
  --mut-scores outputs/NR_000027.1/visualization/mutant_scores.json `
  --base-pairs-wt outputs/NR_000027.1/visualization/wt_base_pairs.json `
  --base-pairs-mut outputs/NR_000027.1/visualization/mut_base_pairs.json `
  --base-pairs-threshold 0.15
```

Utiliser `--base-pairs-mode arcs` pour rendre les connecteurs plus lisibles sur les structures compactes, et `--no-base-pairs` pour desactiver l'affichage depuis un script generique sans supprimer les chemins.

## Backbone discret

Overlay et overlay-delta ajoutent un backbone/tube gris translucide pour suivre la continuite des pseudo-atomes. Il peut etre desactive via `--no-backbone`. Le mode molecule peut afficher le meme tube avec `--backbone`.

## Export PNG

Toutes les commandes acceptent `--export chemin.png` pour capturer une image. Avec `--export-only`, la capture est effectuee puis l'application se ferme sans lancer l'interacteur. Ajouter `--export-scale 2` (ou plus) permet de sur-echantillonner le rendu PNG sans modifier la taille de la fenetre. Les options `--export-view`, `--export-hide-ui` et `--export-background` permettent de figer une camera deterministe et un fond adapte (par exemple `--export-view iso --export-hide-ui --export-background white` pour une figure d'article).

```powershell
python scripts/view_3d.py --mode overlay `
  --wt outputs/runA/visualization/wt_structure.pdb `
  --mut outputs/runA/visualization/mutant_structure.pdb `
  --scores outputs/runA/visualization/mutant_scores.json `
  --export outputs/runA/visualization/overlay_iso.png `
  --export-view iso `
  --export-hide-ui `
  --export-background white `
  --export-only
```
