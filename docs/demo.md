# Démo hors ligne : comparer deux ARN

La démo utilise deux séquences synthétiques de 12 nucléotides, conservées dans
[`data/examples/rna_variant.json`](../data/examples/rna_variant.json).
Le mutant porte une substitution C → A à la position 12, comptée depuis 1.
Ces séquences illustrent le pipeline ; elles ne proviennent ni d'un organisme
identifié ni d'un patient.

Après l'[installation](../README.md#installation), lancer depuis la racine
du dépôt, sous Windows, macOS ou Linux :

```bash
python scripts/run_demo.py
```

Pour choisir le dossier de sortie, y compris un chemin contenant des espaces :

```bash
python scripts/run_demo.py --output-dir "outputs/ma demo"
```

L'exécution n'appelle pas le NCBI, ne nécessite ni clé API ni email et n'ouvre
aucune fenêtre. Elle utilise le même `run_and_export` que la CLI principale.
Les séquences et les paramètres sont fixes : les structures et les scores sont
reproductibles, tandis que les horodatages reflètent chaque exécution.
Relancer dans le même dossier remplace les sorties de cette démo.

## Résultats attendus

Pour les inspecter dans une figure interactive, après régénération avec le code
actuel de `main` :

```bash
python scripts/view_secondary.py --manifest outputs_demo/offline_demo/visualization/visualization_manifest.json
```

La figure utilise les paires exportées dans `secondary_structures.json` et les
scores calculés, sans inventer une conformation 3D. Ajouter
`--export outputs_demo/comparison_secondary.png --export-only` pour un PNG sans fenêtre.
Le rendu par défaut montre les tiges et boucles, les bases A/U/G/C et les extrémités
5′/3′. **Vert** : paire conservée ; **bleu pointillé** : paire perdue ; **rouge** :
nouvelle paire. Ajouter `--view arcs` pour retrouver les arcs et les barres de delta.
L'application desktop ouvre cette même vue 2D en premier ; voir le
[guide desktop](desktop.md#parcours-arn).

Avec les paramètres par défaut :

```text
WT      GGGGAAAACCCC
        ((((....))))
Mutant  GGGGAAAACCCA
        .(((....))).
```

| Mesure | Référence (WT) | Mutant |
| --- | ---: | ---: |
| Nombre de paires prédites | 4 | 3 |
| Score de stabilité heuristique | 2,0 | 1,5 |

La paire entre les positions 1 et 12 disparaît. Le rapport classe les deux
positions dont les scores changent le plus :

| Position (depuis 1) | Δ = score mutant − score WT |
| --- | ---: |
| 12 | −0,700 |
| 1 | −0,600 |

Les dix autres positions ont un delta nul. Une substitution peut donc modifier
le score de son partenaire d'appariement, même si sa base reste identique.

## Fichiers produits

```text
outputs_demo/
  offline_demo.md
  offline_demo.json
  offline_demo_hotspots.csv
  offline_demo/
    visualization/
      impact_summary.json
      hotspots.csv
      wt_structure.pdb
      mutant_structure.pdb
      wt_scores.json
      mutant_scores.json
      visualization_manifest.json
      secondary_structures.json
      wt_base_pairs.json
      mut_base_pairs.json
```

Ouvrir `offline_demo.md` pour lire le tableau des positions les plus modifiées.
Le champ `impact_summary` du rapport JSON et le fichier
`visualization/impact_summary.json` contiennent le même résumé, utilisable
dans d'autres analyses. Les hotspots utilisent des positions comptées depuis 1
(voir le [guide de numérotation](coordinates.md) pour la correspondance détaillée).

Le score par position dépend de l'appariement prédit, de la base et de la
distance entre partenaires. Les deltas mesurent des différences de cette
heuristique, pas des probabilités de pathogénicité ni une énergie physique.
Les fichiers PDB sont des projections géométriques simplifiées ; cette démo
porte sur les appariements secondaires et leurs scores, sans validation
expérimentale ou interprétation clinique.

## Régénérer l'illustration

```bash
python scripts/plot_demo.py
```

Cette commande recalcule les résultats avec le pipeline et génère
`docs/assets/demo_delta.png` avec Matplotlib, sans ouvrir de fenêtre.
Les arcs représentent les paires prédites ; les barres montrent les deltas de
score calculés. L'option `--output "chemin/figure.png"` permet de choisir le
fichier de destination. Matplotlib est déclaré parmi les dépendances du projet.
