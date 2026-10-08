# Benchmark manuel de l'alignement sous Windows

Mesures du 8 octobre 2026 : Intel Core i7-13620H, Windows build 26300, Python 3.11.9, Biopython 1.88.

Une mesure par cas, dans un processus neuf. Seed 42, 10 substitutions, 3 insertions et 3 deletions generees. Aucun repliement ARN. Temps du moteur et detection du second optimum ; pic memoire du processus entier en MiB, imports inclus. Ce ne sont ni des moyennes ni des garanties de performance.

| Profil | Motif | WT / MUT (nt) | Cellules | Temps (s) | Pic (MiB) | Identite | S / I / D trouves | Ambigu |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| 1 M | deterministic-random | 1000 / 1000 | 1,000,000 | 0.0116 | 42.28 | 98.405% | 10 / 3 / 3 | non |
| 1 M | repetitive | 1000 / 1000 | 1,000,000 | 0.0087 | 42.08 | 98.504% | 9 / 3 / 3 | oui |
| 4 M | deterministic-random | 2000 / 2000 | 4,000,000 | 0.0478 | 48.41 | 99.201% | 10 / 3 / 3 | oui |
| 4 M | repetitive | 2000 / 2000 | 4,000,000 | 0.0308 | 48.26 | 99.201% | 10 / 3 / 3 | oui |
| 10 M | deterministic-random | 3162 / 3162 | 9,998,244 | 0.1059 | 59.70 | 99.494% | 10 / 3 / 3 | oui |
| 10 M | repetitive | 3162 / 3162 | 9,998,244 | 0.0675 | 59.89 | 99.494% | 10 / 3 / 3 | non |
| 25 M | deterministic-random | 5000 / 5000 | 25,000,000 | 0.3352 | 89.48 | 99.680% | 10 / 3 / 3 | oui |
| 25 M | repetitive | 5000 / 5000 | 25,000,000 | 0.1648 | 89.38 | 99.700% | 9 / 3 / 3 | oui |
| 50 M | deterministic-random | 7070 / 7070 | 49,984,900 | 0.5411 | 137.55 | 99.774% | 10 / 3 / 3 | oui |
| 50 M | repetitive | 7070 / 7070 | 49,984,900 | 0.3267 | 137.47 | 99.774% | 10 / 3 / 3 | oui |

Toutes ces demandes ont ete acceptees avec la limite explicitement indiquee. Les profils 10 M et 50 M utilisent respectivement 3162 et 7070 bases : le produit reste legerement inferieur a la limite.

[Données brutes et configurations](alignment_windows.json). Voir [les commandes et la méthode](../configuration.md).

Le cas aleatoire peut lui aussi avoir des optima ex aequo : repetitions locales et indels suffisent. Le nombre de differences retrouvees ne constitue pas une reconstruction certaine des operations de generation.

Aucun benchmark Linux mesure dans cette passe. La CI teste les bornes sur des cas minuscules ; les grandes matrices restent des diagnostics manuels. Ces dix mesures ne justifient pas de relever le defaut de production : il reste a 4 M.
