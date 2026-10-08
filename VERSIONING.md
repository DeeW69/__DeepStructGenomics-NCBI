# Politique de versionnement

La version de développement est **0.5.0rc1**, définie dans [pyproject.toml](pyproject.toml),
affichée **v0.5.0-rc1** dans la documentation. La dernière release stable reste **v0.4.6**.
Le [changelog](CHANGELOG.md) décrit les changements et les
[releases GitHub](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases)
regroupent les publications et leurs artefacts.

## Numérotation

Le format est `MAJOR.MINOR.PATCH`, précédé de `v` pour les tags Git.

| Niveau | Usage habituel |
| --- | --- |
| PATCH | Corrections, stabilisation et documentation sans rupture. |
| MINOR | Fonctionnalités et formats supplémentaires. |
| MAJOR | Nouvelle génération avec changements incompatibles majeurs. |

Le projet reste en série `0.x`. Toute modification incompatible doit être annoncée
et documentée ; elle ne doit pas être introduite silencieusement dans un correctif.
`v0.4.x` désigne une série, **jamais un numéro de paquet ou de tag**.

### Exception convenue pour la série v0.4.x

Les six points de développement utilisent **v0.4.1 à v0.4.6**, avec un commit
fonctionnel par point, même lorsqu'ils ajoutent des fonctionnalités.

| Commits / incréments | Publication globale |
| --- | --- |
| v0.4.1 : cache ; v0.4.2 : lots ; v0.4.3 : évaluation ; v0.4.4 : synthèse Hi-C | v0.4.4 |
| v0.4.5 : benchmark ARN expérimental ; v0.4.6 : Hi-C avancé | v0.4.6 |

Les incréments intermédiaires ne nécessitent pas de releases GitHub séparées.
La publication v0.4.6 inclut tous les changements antérieurs sans déplacer le
tag v0.4.4. Hors de cette série, les règles habituelles s'appliquent.

## Éléments à synchroniser

| Élément | Publication actuelle |
| --- | --- |
| Version Python | `0.4.6` |
| Tag Git annoté | `v0.4.6` |
| Titre de release | `DeepStructGenomics-NCBI v0.4.6` |
| Notes | `docs/releases/v0.4.6.md` |

Avant publication, aligner le paquet, le README, le changelog, ce document et
les notes de release. Les changements postérieurs sont consignés dans une section
« Non publié » jusqu'à la prochaine version. Un numéro de paquet sur `main` ne
suffit pas à identifier une révision : utiliser le tag ou le SHA du commit.

Si les dépendances changent, mettre à jour leurs déclarations et les verrous
concernés dans le même commit. Le socle utilise `requirements-lock.txt` ;
ViennaRNA/SciPy utilisent `requirements-research-lock.txt`, contraint par le socle.
L'interface desktop facultative utilise `requirements-gui-lock.txt`. Sa première
version est préparée pour v0.5.0-rc1 ; les changements sur `main` restent « Non publié ».
Voir le [guide des dépendances](docs/dependencies.md).

## Validation et publication

1. Installer les verrous dans un environnement neuf ; exécuter `pip check`, les
   tests, les démos, le benchmark ARN et l'analyse Hi-C avec export PNG.
2. Commiter les changements. Construire wheel et archive source depuis ce commit,
   puis vérifier le paquet installé hors du chemin source.
3. Publier les commits et attendre la réussite des jobs GitHub Actions Linux /
   Python 3.10 et Windows / Python 3.12, avec et sans extensions de recherche.
4. Créer un tag annoté sur **ce commit validé**, puis publier la release associée.
5. Joindre wheel, archive source, archive de démonstration/benchmark et
   `SHA256SUMS.txt`. Vérifier les empreintes des fichiers téléchargés après publication.

Les archives sources incluent les verrous, scripts, corpus réduit et documentation.
Les résultats scientifiques publiés indiquent données, provenance, versions,
paramètres, limites et conditions de reproduction. Une amélioration logicielle
ne vaut pas validation biologique ; les méthodes exploratoires sont identifiées.
Le workflow CI ne publie automatiquement ni release ni paquet PyPI.

La [checklist](docs/release_checklist.md) complète cette procédure.

## Préversions et corrections

Pour cette candidate : paquet PEP 440 `0.5.0rc1`, futur tag `v0.5.0-rc1`, release
marquée **Pre-release** et limites explicites. Ni tag ni release ne sont créés
automatiquement par la préparation du paquet. Voir les [notes de préparation](docs/releases/v0.5.0-rc1.md).

Une version publiée conserve son tag et ses artefacts. Toute correction ultérieure
utilise un nouveau numéro ; ne pas déplacer un tag publié ni réutiliser une version
retirée. Les évolutions CLI/API et les conventions de coordonnées sont documentées
dans les notes, avec une procédure de migration en cas de rupture.

## Publications

| Version | Date | Contenu |
| --- | --- | --- |
| [v0.4.6](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.4.6) | 7 octobre 2026 | Benchmark ARN expérimental, Hi-C avancé et documentation clarifiée. |
| [v0.4.4](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.4.4) | 7 octobre 2026 | Cache, lots FASTA, évaluation ARN et synthèse Hi-C. |
| [v0.3.0](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.3.0) | 6 octobre 2026 | Ouverture par manifeste, CSV et seuils configurables. |
| [v0.2.1](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.2.1) | 6 octobre 2026 | Dépendances verrouillées, robustesse NCBI et coordonnées. |
| [v0.2.0](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.2.0) | 6 octobre 2026 | FASTA, rapports enrichis, démo et viewers. |

Mainteneur : [DeeW69](https://github.com/DeeW69).
