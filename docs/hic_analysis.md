# Analyse Hi-C avancée

Depuis v0.4.6, `scripts/explore_hic.py --advanced` combine un équilibrage des
contacts, des tests de boucles/frontières et une reconstruction exploratoire.
Les méthodes sont des baselines de recherche explicites, pas des réimplémentations
validées de HiCCUPS, Fit-Hi-C ou d'un pipeline complet de TADs.

## Lancer l'exemple

Installer le [socle et les extensions de recherche](dependencies.md), puis :

```bash
python scripts/explore_hic.py --contacts data/examples/hic_region.tsv --assembly synthetic --bin-size 10000 --advanced --chromosome chrSynthetic --start 0 --end 320000 --missing-as-zero --output-dir outputs_hic_advanced
```

Le dossier de sortie doit être nouveau pour ne pas conserver de tables anciennes
après une analyse en échec. `--no-plot` désactive la figure PNG.
Le mode descriptif de v0.4.4, sans `--advanced`, reste disponible sans SciPy.

## Entrées et coordonnées

TSV UTF-8 : `chrom1`, `start1`, `chrom2`, `start2`, `count`, séparés par tabulations.
Comptes bruts entiers, non négatifs, au plus 10¹². Chaque paire apparaît une seule
fois ; les doublons, y compris leur symétrique, sont refusés. Les coordonnées sont
depuis 0, alignées sur la résolution, avec intervalles `[start, start + bin_size)`.
La région cis est explicite et limitée à **400 bins** pour les calculs denses.
Les lignes hors région et trans sont comptées puis exclues ; leur nombre est rapporté.

`--missing-as-zero` déclare que les contacts absents sont des **zéros observés**.
Ne pas l'utiliser avec un fichier ne contenant que des contacts présélectionnés
ou significatifs : les tests exigent le fond complet. Les bins sans couverture
hors diagonale sont masqués et ne deviennent pas des positions 3D artificielles.
Le programme ne vérifie pas les longueurs chromosomiques de l'assemblage déclaré.

## Équilibrage

Pour une matrice brute C, le programme cherche des poids w tels que les sommes
des lignes de `B[i,j] = w[i] C[i,j] w[j]` valent 1. Les diagonales sont exclues.
La mise à jour multiplicative est amortie par une racine carrée ; les poids,
bins masqués, itérations et erreur relative maximale sont exportés.

Valeurs par défaut : `--tolerance 1e-6`, `--max-iterations 2000`.
Une matrice sans support suffisant peut ne pas converger. Dans ce cas, le JSON
et le rapport diagnostique sont écrits, **aucun test ni modèle 3D n'est calculé**,
et la CLI retourne 2. Les poids sont propres à cette région, pas au chromosome entier.

Ce principe de visibilité égale est celui utilisé par les approches de type
[ICE](https://hiclib.github.io/iced/modules/generated/iced.normalization.ICE_normalization.html).
Notre solveur amorti est distinct du paquet `iced` ; aucune équivalence numérique
avec son implémentation ni correction de tous les biais biologiques n'est revendiquée.

## Boucles candidates

Les comptes bruts sont testés avec une loi binomiale négative, à variance
`mu + alpha * mu²` ; une dispersion nulle utilise la limite de Poisson.
La [paramétrisation SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.nbinom.html)
sert au calcul de la probabilité de dépassement.

Pour chaque pixel, les autres pixels de **même distance génomique**, zéros inclus,
forment le fond. L'exposition est `1 / (w[i] w[j])`. Le taux est la somme des comptes
du fond divisée par la somme des expositions ; la dispersion est estimée par moments
sur les taux corrigés, après soustraction de la variance d'échantillonnage Poisson.
Le pixel testé est exclu de l'estimation du taux et de la dispersion, mais participe
à l'équilibrage global. Au moins 5 autres pixels et 20 comptes de fond sont requis.

La correction Benjamini-Hochberg porte sur **tous les pixels testés de la région**.
Les candidats ont q ≤ `--fdr` (0,05 par défaut) et enrichissement > 1. Les distances
inférieures à `--min-separation` (2 bins) sont ignorées ; les strates trop petites
sont signalées comme non testées. La table conserve aussi les tests non significatifs.

Les paramètres de fond sont estimés, les contacts sont dépendants et aucun réplicat
n'est utilisé : ces q-values sont **exploratoires**, sans garantie de calibration
biologique du taux de faux positifs. Ce module ne remplace pas un caller validé.

## Frontières et domaines candidats

Un score d'isolation mesure la moyenne des contacts dans un carré traversant
chaque frontière, avec `--window-bins 3` de chaque côté. Seules les fenêtres complètes
et non masquées sont testées. Le score brut et son log2 relatif à la médiane positive
sont conservés ; les logarithmes indéfinis valent `null`.

Le modèle nul permute les valeurs entre pixels de même distance génomique,
zéros inclus : `--permutations 999`, `--seed 46`. La p-value basse est
`(1 + nombre de scores permutés <= score observé) / (1 + permutations)`.
Les q-values BH sont calculées sur toutes les frontières testées, séparément des
boucles. Parmi les minima locaux significatifs, le meilleur est retenu par fenêtre.
Les intervalles entre ces frontières constituent des **domaines candidats**.

Les permutations préservent les distributions par distance, mais ni les marges
équilibrées ni la dépendance spatiale. Le nombre de permutations limite la résolution
des p-values. Des domaines sans frontière significative ne prouvent pas son absence.

## Reconstruction 3D

Les contacts équilibrés positifs deviennent des distances `contact^(-exposant)`
avec `--distance-exponent 0.3333333333333333`. Les
[plus courts chemins](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.shortest_path.html)
complètent les distances dans le graphe connecté. Une mise à l'échelle
multidimensionnelle classique (MDS) projette ces distances en trois dimensions.
Le programme ne relie pas artificiellement les composantes déconnectées :
`status=disconnected`, aucune coordonnée, les autres résultats restent disponibles.

Les unités sont arbitraires (médiane des distances cibles = 1). Le stress est
`sqrt(sum((distance_3D - distance_cible)²) / sum(distance_cible²))` ; la fraction
de valeurs propres négatives indique le caractère non euclidien des distances.
Ces diagnostics mesurent l'ajustement aux hypothèses, pas la fidélité biologique.
Rotation/réflexion des coordonnées et petites différences numériques entre systèmes
sont possibles. Il ne s'agit pas d'une structure atomique ni d'une mesure de distances.

## Sorties

| Fichier | Contenu |
| --- | --- |
| `hic_analysis.json` | Provenance SHA-256, versions, paramètres, diagnostics et résultats complets. |
| `hic_analysis.md` | Synthèse et limites de l'analyse. |
| `balanced_contacts.tsv` | Poids équilibrés hors diagonale, triangle supérieur. |
| `loop_tests.tsv` | Comptes, fond, dispersion, p/q et sélection. |
| `boundary_tests.tsv` | Scores d'isolation, p/q et frontières retenues. |
| `candidate_domains.bed` | Intervalles candidats depuis 0, borne droite exclue. |
| `coordinates.tsv` | Bins actifs et coordonnées inférées en unités arbitraires. |
| `hic_analysis.png` | Carte des contacts équilibrés et projection 3D. |

## Validation livrée

L'exemple synthétique contient 32 bins, deux domaines séparés au bin 16 et une
boucle ajoutée entre les bins 5 et 10. Le pipeline retrouve cette frontière et
cette seule boucle avec les paramètres par défaut. Le générateur et sa vérité
de simulation sont inclus ; ce cas ne constitue pas une validation Hi-C expérimentale.
Le seuil et le modèle ne doivent pas être optimisés sur ce seul exemple.

Les tests couvrent aussi biais connus, bins vides, support non convergent,
probabilités BH connues, reproductibilité des permutations et reconstruction
d'une géométrie euclidienne connue à rotation/échelle près.
