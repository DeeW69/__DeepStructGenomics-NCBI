# Alignement WT/MUT — v0.5.0-rc1

L'alignement global est activé par défaut lorsqu'un mutant est fourni. Il met en
correspondance les séquences avant la comparaison des résultats de repliement.
Il ne prédit ni structure secondaire, ni conformation 3D, ni effet biologique.

## Modèle minimal

`alignment/pairwise.py` appelle `Bio.Align.PairwiseAligner` dans le processus de
calcul. Biopython est déjà une dépendance du socle ; aucun service externe n'est
appelé. Les valeurs par défaut sont centralisées dans `AlignmentConfig` et
surchargeables par [.env, CLI et GUI](configuration.md) :

| Paramètre | Valeur |
| --- | --- |
| Match | +2 |
| Substitution | −1 |
| Ouverture de gap | −5 |
| Extension de gap | −1 |

Les gaps terminaux ont les mêmes pénalités que les gaps internes. Un gap de
longueur *n* coûte `−5 − (n − 1)`, selon la définition des
[scores affines de Biopython](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html).
Ces valeurs constituent une convention pour des variantes proches, pas des
paramètres biologiquement universels. Une substitution peut être préférée à deux
indels ; le logiciel ne prétend pas retrouver l'histoire réelle d'une mutation.

Le moteur prend le premier traceback optimal sur les séquences renversées puis
remet les colonnes dans le sens original. Cela place à droite le gap des répétitions
simples : `ACGU → ACGGU` donne `ACG-U / ACGGU` (insertion G4), et l'inverse une
délétion G4. Il vérifie seulement l'existence d'un second optimum, sans énumérer
tous les alignements. Les ex æquo donnent un avertissement. L'emplacement d'un
indel dans une répétition reste ambigu ; il ne faut pas interpréter cette règle
comme une preuve d'homologie ou une normalisation HGVS.

L'identité est **matches / nombre total de colonnes, gaps inclus**. Les nombres
de substitutions, insertions et délétions comptent des bases, pas des événements.
Le CIGAR utilise `=`, `X`, `I`, `D` depuis la référence vers le mutant.

## Coordonnées et source de vérité

`AlignmentResult` contient les deux chaînes avec gaps et les colonnes dérivées.
Les maps réciproques et les maps vers les colonnes sont calculées depuis ces
colonnes, jamais par une seconde heuristique dans Qt.

| Champ | Convention |
| --- | --- |
| `AlignmentColumn.index` | Colonne depuis 1 |
| `reference_position`, `mutant_position` | Position dans sa séquence depuis 1, ou `null` pour un gap |
| `ref_to_mut`, `mut_to_ref` | Position depuis 1 → position depuis 1 ou `None` |
| `ref_to_alignment`, `mut_to_alignment` | Position depuis 1 → colonne depuis 1 |
| Paires dans `structure.base_pairs` et `secondary_structures.json` | Indices natifs depuis 0, inchangés |
| Paires affichées sur les structures 2D | Positions natives depuis 1 |
| Axe du delta et liste des changements de paires | Colonnes depuis 1 si alignement actif |

Un clic sur WT ou MUT sélectionne la colonne puis l'homologue s'il existe. Un
gap ne reçoit ni base fictive, ni score zéro, ni delta. L'inspecteur affiche les
deux positions natives et la colonne ; la mini-séquence comprend jusqu'à quatre
colonnes de chaque côté. La table **Alignement WT/MUT** donne les positions, bases,
variations, scores et contextes, avec défilement horizontal.

La géométrie 2D reste propre aux structures lorsque des gaps existent. Les couleurs
des paires et des deltas utilisent les correspondances. En 3D schématique, une base
mutante sans homologue est grise ; son clic ouvre la bonne colonne. Les connecteurs
entre structures (`--links`) sont refusés avec un alignement : leur ancien calcul
de distance par indices bruts n'est pas adapté aux indels. Les coordonnées 3D ne
sont jamais réinterprétées comme une déformation moléculaire.
Les options avancées de couleur/seuil des paires 3D restent aussi indisponibles
avec des indels ; le mode lignes de la GUI est conservé.

## Comparaison des paires

Une paire est conservée si ses deux extrémités correspondent à une paire mutante.
Les catégories sont disjointes :

- **Perdue** : deux homologues existent, mais la paire n'existe plus.
- **Nouvelle** : deux homologues existent, mais la paire n'existait pas dans WT.
- **Supprimée avec une base** : la paire WT implique au moins une base supprimée.
- **Nouvelle avec insertion** : la paire mutante implique au moins une base insérée.

Un changement de partenaire peut donc afficher perte et gain ensemble. Les cartes
de synthèse comptent toutes les paires absentes/nouvelles, indels inclus, avec un
libellé explicite lorsqu'un gap intervient. Les catégories sont détaillées dans
**Delta structural** et l'inspecteur. La légende des dessins regroupe en bleu
les paires absentes du mutant et en rouge les nouvelles, y compris celles liées
aux gaps ; le libellé « perdue/supprimée » l'explicite.

## Limites et fallback explicite

Par défaut, avant tout repliement, la demande est refusée si une séquence dépasse
10 000 bases ou si `longueur WT × longueur MUT > 4 000 000`. Les limites sont
configurables. Ce nombre théorique de cellules est un indicateur de travail,
pas une mesure de RAM allouée. Il n'accélère pas Nussinov : les longs ARN peuvent rester
coûteux à replier. Le worker reste annulable et la progression est indéterminée.

Il n'y a pas de fallback silencieux. Pour une comparaison brute, décocher
**Activer l'alignement WT/MUT** dans les paramètres avancés, ou passer `--no-alignment`
en CLI. La politique `--overflow-policy positional` autorise explicitement un repli
par position en cas de dépassement ; `reject` reste le défaut.
La GUI affiche **Alignement non calculé — comparaison brute par position, sans alignement**.
Les anciens résultats sans alignement utilisent ce même mode, sans modification
des fichiers et sans recalcul en thread GUI. Un alignement enregistré incohérent
est refusé avec une erreur récupérable ; il n'est jamais remplacé silencieusement.

## Exports et migration

Le bloc `alignment` est inclus dans le rapport JSON et dans
`secondary_structures.json` référencé par le manifeste : paramètres, implémentation
et version Biopython, règle de départage, score, identité, counts, CIGAR, chaînes,
colonnes, booléen `ambiguous` et avertissements. Le bloc `alignment_run` conserve
la configuration effective, les longueurs natives, les cellules et le statut,
y compris lorsque l'alignement est désactivé ou que le repli positionnel est choisi.
Aucun fichier séparé n'est nécessaire pour rouvrir
l'analyse. Les maps sont reconstruites depuis ce bloc et vérifiées à la lecture.

Les clés historiques restent présentes, mais une analyse **alignée** déclare :

- `variant_analysis.coordinate_system = "alignment_column_0based"` : le champ
  historique `substitutions[].position` est l'indice de colonne depuis 0. Cette
  liste historique contient toutes les différences, indels inclus ; consulter
  `alignment.columns` pour le type et les positions natives depuis 1.
- `impact_summary.coordinate_system = "alignment_column_1based"` : `hotspots[].position`
  et la colonne `position` des CSV hotspots désignent une colonne d'alignement.
  Les bases exportées sont celles de cette colonne et un gap n'a pas de hotspot.

Les rapports Markdown l'indiquent ; les chaînes alignées sont affichées seulement
jusqu'à 120 colonnes (les données complètes restent dans le JSON). Pour
retrouver une position WT/MUT, utiliser `alignment.columns[position - 1]` pour un
hotspot. Les scores individuels, fichiers PDB et fichiers de paires conservent
leurs coordonnées natives. Les consommateurs des anciens rapports doivent lire
`coordinate_system` avant de traiter de nouvelles analyses ; son absence conserve
le sens positionnel historique. `--no-alignment` garde ce mode de comparaison.

## Exemple reproductible

```powershell
python scripts/run_pipeline.py --sequence GGGGAAAACCCC --mutant-sequence GGGGUAAAACCCA --label demo_indel --output-dir outputs_indel
python scripts/run_gui.py --manifest outputs_indel/demo_indel/visualization/visualization_manifest.json
```

La démo contient une insertion U5 et la substitution WT C12 → MUT A13. Les positions
suivantes restent mises en correspondance. L'algorithme peut signaler une ambiguïté
sur d'autres entrées, en particulier répétitives ou éloignées. Les profils GUI et
les commandes de benchmark sont décrits dans [configuration.md](configuration.md).
