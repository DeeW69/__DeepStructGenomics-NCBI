# Politique de versionnement

Ce document définit les règles de versionnement et de publication de
DeepStructGenomics-NCBI. Les changements livrés sont décrits dans le
[changelog](CHANGELOG.md) et les [releases GitHub](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases).

## Numérotation

Le projet utilise trois nombres : `MAJOR.MINOR.PATCH`.

| Niveau | Utilisation dans ce projet | Exemple |
| --- | --- | --- |
| PATCH | Correction, stabilisation ou documentation sans rupture des interfaces existantes. | `0.2.0` → `0.2.1` |
| MINOR | Fonctionnalité ou format de sortie supplémentaire. | `0.2.1` → `0.3.0` |
| MAJOR | Nouvelle génération avec changements incompatibles majeurs. | `1.0.0` → `2.0.0` |

La série `0.x` correspond à un projet encore en développement. Les évolutions
incompatibles pendant cette phase doivent être annoncées dans une version
mineure et expliquées dans les notes de migration ; elles ne doivent pas être
introduites silencieusement dans une version corrective.

Exception explicite pour la série v0.4.x : les quatre points de la roadmap
utilisent les incréments v0.4.1 à v0.4.4, un commit par point. La release globale
est publiée sous le tag précis `v0.4.4` et le titre « v0.4.x — exploration et lots ».
Les incréments intermédiaires ne nécessitent pas de releases GitHub séparées.

`v0.2.x` désigne une **série**, pas une version publiable. Une release utilise
toujours un numéro précis, par exemple `v0.2.1`.

## Version du paquet, tag et release

La valeur `[project].version` de [pyproject.toml](pyproject.toml) est la référence
pour la version embarquée dans le wheel et l'archive source.

| Élément | Exemple pour la stabilisation publiée |
| --- | --- |
| Version Python dans `pyproject.toml` | `0.2.1` |
| Tag Git annoté | `v0.2.1` |
| Titre de release | `DeepStructGenomics-NCBI v0.2.1` |
| Notes de release | `docs/releases/v0.2.1.md` |

Les tags identifient les commits publiés. La branche `main` peut contenir des
changements ultérieurs : consulter le tag pour reproduire exactement une
release. Le numéro du paquet est mis à jour lors de la préparation de la
prochaine publication ; il ne suffit donc pas à identifier chaque commit de `main`.

Pour une préversion, utiliser une version Python telle que `0.3.0a1`,
`0.3.0b1` ou `0.3.0rc1`, avec un tag correspondant préfixé par `v`.
Cocher **Pre-release** sur GitHub et préciser les limites de validation.

## Fichiers à synchroniser

Avant chaque publication :

1. Mettre à jour `[project].version` dans [pyproject.toml](pyproject.toml).
2. Regrouper les changements dans une entrée datée du [changelog](CHANGELOG.md).
3. Rédiger les notes dans `docs/releases/vX.Y.Z.md` : changements, compatibilité,
   installation, validation et fichiers disponibles.
4. Aligner le badge, la version publiée, les features et la roadmap du
   [README](README.md).
5. Si les dépendances changent, synchroniser leurs déclarations et régénérer
   le verrou selon le [guide des dépendances](docs/dependencies.md).

Les changements ajoutés après une release sont consignés dans une section
**Non publié** du changelog jusqu'à la préparation de la version suivante.
Une simple modification de documentation n'impose pas une release immédiate.

## Validation et publication

La [checklist de release](docs/release_checklist.md) détaille les commandes.
L'ordre de publication est le suivant :

1. Valider dans un environnement neuf installé depuis `requirements-lock.txt` :
   compatibilité avec `pip check`, tests, démo hors ligne et export PNG.
2. Commiter les changements de version et de documentation. Construire le wheel
   et l'archive source depuis ce commit, puis vérifier le paquet installé.
   Si une correction est nécessaire, créer un nouveau commit et reconstruire
   les artefacts concernés avant de poursuivre.
3. Pousser le commit et attendre la réussite des deux jobs GitHub Actions :
   Linux / Python 3.10 et Windows / Python 3.12.
4. Créer un tag annoté sur **ce même commit validé**, puis pousser ce tag.
5. Préparer une release GitHub associée au tag, reprendre les notes et joindre
   les fichiers vérifiés avant de publier.

Exemple de commandes à adapter à une **nouvelle** version et à son commit :

```bash
git tag -a vX.Y.Z COMMIT_VALIDE -m "Release vX.Y.Z"
git push origin vX.Y.Z
```

Fichiers à joindre : wheel, archive source, archive de démonstration et
`SHA256SUMS.txt` contenant leurs empreintes SHA-256. L'archive source inclut
le verrou et les guides. Installer uniquement le wheel n'applique pas
automatiquement le verrou des dépendances.

Le workflow actuel exécute les tests ; il ne publie pas automatiquement une
release GitHub ni un paquet sur PyPI.

## Compatibilité et corrections après publication

Les modifications de la CLI, de l'API Python ou des formats de rapport doivent
être décrites explicitement. En particulier, les champs JSON historiques de
positions conservent leurs conventions ; voir le
[guide des coordonnées](docs/coordinates.md).

Une version publiée conserve son tag et ses artefacts. Pour corriger une
release, publier une nouvelle version corrective et signaler le problème dans
ses notes. Une suppression exceptionnelle est décrite dans la checklist ;
ne pas réutiliser le numéro retiré ni déplacer un tag déjà publié.

## Historique des releases

| Version | Date | Contenu |
| --- | --- | --- |
| [v0.4.4](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.4.4) | 7 octobre 2026 | Série v0.4.x : cache, lots FASTA, évaluation ARN et exploration Hi-C. |
| [v0.2.0](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.2.0) | 6 octobre 2026 | Entrées FASTA, rapports enrichis, démo et améliorations du viewer. |
| [v0.2.1](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.2.1) | 6 octobre 2026 | Dépendances verrouillées, robustesse NCBI, 132 tests et guide des coordonnées. |
| [v0.3.0](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.3.0) | 6 octobre 2026 | Ouverture par manifeste, exports CSV et seuils configurables. |

Mainteneur : [DeeW69](https://github.com/DeeW69).
