# Historique des versions

## Non publié — préparation v0.5.0

- Vue ARN principale en tiges/boucles 2D, bases colorées, extrémités 5′/3′ et
  positions communes lorsque les appariements sont compatibles. Comparaison
  des paires conservées (vert), perdues (bleu) et nouvelles (rouge).
- Onglets Structure secondaire, Séquence, Delta structural et 3D schématique ;
  inspection du contexte structural, sélection synchronisée et lettres adaptées au zoom.
- Vue VTK complémentaire avec glyphes plus petits/mats et lettres des bases.
- CLI `view_secondary.py` : tiges/boucles par défaut, `--view arcs` pour la vue précédente.

- Interface desktop PySide6 facultative : accueil, séquences, comparaison et Hi-C.
- Réutilisation du pipeline ARN et de l'analyse Hi-C dans un processus séparé,
  annulable, avec un dossier unique par calcul et un historique local des résultats.
- Comparaison 2D et scores par position, inspection des bases/paires, scène VTK
  intégrée à la demande ; informations techniques repliables et géométrie explicitée.
- Heatmap Hi-C, boucles et reconstruction inférée depuis les rapports du moteur.
- Figures PNG, réouverture des manifests/rapports et accès aux exports existants.
- Extra `gui`, verrou à empreintes séparé, lanceur installé `deepstructgenomics`
  et tests des parcours desktop hors écran sous Linux/Windows.


## 0.4.6 — 7 octobre 2026 — release globale

- Analyse Hi-C cis : équilibrage à visibilité égale, tests binomiaux négatifs de
  boucles, permutations d'isolation pour les frontières et correction BH.
- Reconstruction par distances de contacts et MDS, avec stress, masques et
  diagnostic des graphes déconnectés ; sorties JSON/Markdown/TSV/BED/PNG.
- Exemple synthétique reproductible et tests numériques, dont non-convergence
  et surdispersion. Limite de 400 bins et hypothèses documentées.
- README nettoyé, guides pipeline/Hi-C séparés, VERSIONING.md actualisé.
- Inclut v0.4.5 ; 197 tests locaux avec extensions de recherche.

## 0.4.5 — 7 octobre 2026

- Benchmark sur 30 références bpRNA/PDB, provenance et sélection reproductible.
- Exécution optionnelle de ViennaRNA 2.7.2 et archivage des prédictions/métriques.
- Dépendances de recherche optionnelles avec verrou séparé.

## 0.4.4 — 7 octobre 2026 — release globale v0.4.x

- Import exploratoire Hi-C TSV avec assemblage/résolution explicites, contrôles
  des coordonnées et des doublons ; rapports JSON/Markdown, poids cis/trans,
  couverture par bin, répartition des distances et contacts de poids maximal.
- Inclut les incréments v0.4.1, v0.4.2 et v0.4.3 ci-dessous.
- Ni correction des biais Hi-C ni reconstruction 3D ; exemples synthétiques.

## 0.4.3 — 7 octobre 2026

- Évaluation ARN : précision/rappel/F1, méthodes témoins et prédictions externes.
- Jeu synthétique reproductible ; validation biologique indépendante encore nécessaire.

## 0.4.2 — 7 octobre 2026

- FASTA par lots : rapports isolés, erreurs par entrée et synthèse JSON.

## 0.4.1 — 7 octobre 2026

- Cache local NCBI atomique, provenance SHA-256 et rafraîchissement explicite.
- CLI : `--cache-dir` et `--refresh-cache`.

### Visualisation incluse depuis v0.3.0

- Vue interactive et export PNG des appariements secondaires via
  `scripts/view_secondary.py --manifest`, à partir des prédictions exportées :
  référence/mutant côte à côte, paires perdues/gagnées, bases différentes et deltas.
- Nouvel artefact `secondary_structures.json` et exports de paires depuis 1,
  détectés automatiquement par le viewer VTK. Régénérer les anciens résultats
  pour utiliser la nouvelle vue secondaire.
- Correction de la légende VTK surdimensionnée, palette à zéro neutre clair et
  état des connecteurs synchronisé avec leur affichage.
- Documentation corrigée : la géométrie pseudo-3D dépend seulement de la
  longueur, elle ne prédit pas une conformation moléculaire.
- Dépendance directe Matplotlib déclarée, sans mise à jour des versions verrouillées.

## 0.3.0 — 6 octobre 2026

[Release GitHub et artefacts](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.3.0)

- Ouverture simplifiée d'un résultat 3D avec l'option `--manifest` dans `scripts/view_3d.py` :
  détection automatique du viewer (`overlay-delta`, `overlay`, `molecule`) et résolution
  automatique des chemins de structures et de scores relatifs au dossier du manifeste.
- Export tabulaire des positions les plus modifiées (hotspots) au format CSV
  (`<run_id>_hotspots.csv` dans le dossier des rapports et `hotspots.csv` dans
  les artefacts de visualisation) avec colonnes `position`, `reference`, `mutant`,
  `delta` et `abs_delta` pour l'analyse dans des tableurs ou notebooks.
- Nouvelles options CLI dans `scripts/run_pipeline.py` :
  `--top-k` pour contrôler le nombre maximal de hotspots affichés,
  `--min-abs-delta` (alias `--delta-threshold`) pour fixer le seuil minimal de variation,
  et `--base-pair-threshold` pour ajuster la sensibilité du comptage des paires affectées.
- Ajout de `VERSIONING.md` : règles de numérotation, synchronisation des fichiers,
  préversions et publication des releases. Ce guide est lié depuis le README
  et inclus dans les futures archives source.

## 0.2.1 — 6 octobre 2026

[Release GitHub et artefacts](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.2.1)

- Tests hors réseau du client NCBI : réponses FASTA vides ou invalides,
  résumés de métadonnées vides ou mal formés, erreurs HTTP 400/404/429/500 et
  délais dépassés sur EFetch et ESummary.
- Correction de la lecture des résumés JSON contenant `null`, une liste ou une
  entrée de métadonnées non objet : la séquence valide est conservée avec des
  métadonnées `{}`, comme pour une réponse vide ou un JSON invalide.
- Tests de la CLI pour les erreurs HTTP, de connexion et de délai : arrêt avant
  la prédiction et l'export, message lisible sans afficher l'email ni la clé API.
  Les erreurs HTTP et réseau continuent à se propager dans l'API
  Python ; elles ne sont pas traitées comme des métadonnées absentes.
- Verrou universel `requirements-lock.txt` pour les dépendances du pipeline et
  des tests, avec versions exactes, marqueurs Python/plateforme et empreintes
  SHA-256 des distributions.
- Installation à partir du verrou dans les jobs Linux / Python 3.10 et Windows /
  Python 3.12 ; contrôle des empreintes, roues précompilées et `pip check`.
- Collecte pytest limitée à `tests/` pour éviter les doublons si des copies du
  code sont présentes dans les dossiers de sortie ou de construction.
- Guide d'installation et de régénération dans `docs/dependencies.md` ; le verrou
  est inclus dans les futures archives source. Les plages de dépendances du
  paquet restent définies dans `pyproject.toml`.
- Documentation de la numérotation des positions dans
  `docs/coordinates.md` : correspondance entre conventions 0-based
  (algorithmes Python, champs historiques JSON) et 1-based (validation des
  séquences, rapports Markdown, hotspots, résidus PDB, scores et viewer VTK),
  avec tableau de synthèse, exemple pas à pas sur la démo 12 nt et formules de conversion.

## 0.2.0 — 6 octobre 2026

[Release GitHub et artefacts](https://github.com/DeeW69/__DeepStructGenomics-NCBI/releases/tag/v0.2.0)

### Entrées FASTA et validation

- Options `--fasta` et `--mutant-fasta` pour charger une référence et un mutant
  depuis des fichiers locaux contenant chacun une seule séquence.
- Lecture des séquences réparties sur plusieurs lignes, des minuscules et des
  espaces ; conservation de l'identifiant, de la description et de la provenance
  FASTA. L'option `--label` permet de nommer les sorties de la référence.
- Exemples FASTA synthétiques fournis dans `data/examples/`, utilisables hors
  ligne avec la CLI principale.
- Normalisation ADN → ARN (`T` → `U`) appliquée aux séquences NCBI, directes et
  FASTA avant les calculs et la comparaison référence/mutant.
- **Compatibilité :** les caractères hors `A`, `C`, `G`, `T`, `U` (codes ambigus,
  gaps, chiffres, etc.) déclenchent désormais une erreur indiquant leur position,
  comptée depuis 1 sans les espaces. Ils ne sont plus supprimés silencieusement,
  ce qui évite de décaler les positions analysées. Les séquences vides et les
  fichiers FASTA mal formés ou contenant plusieurs séquences sont refusés.

### Comparaison référence/mutant

- Résumé des variations par position dans le rapport Markdown : statistiques
  des deltas, jusqu'à dix positions les plus modifiées et comptage des paires
  contenant une position au-dessus du seuil.
- Un même résumé est disponible dans `PipelineResult.impact_summary`, dans le
  champ `impact_summary` du rapport JSON et dans `visualization/impact_summary.json`.
- L'absence de mutant est signalée explicitement et distinguée d'une comparaison
  dont les scores sont identiques. Les comparaisons de longueurs différentes
  précisent qu'aucun alignement n'est effectué.
- Démo synthétique hors ligne : `python scripts/run_demo.py`, avec résultats
  attendus et illustration dans le README.

### Visualisation

- Exports PNG avec presets de caméra, choix du fond, masquage des éléments
  d'interface et facteur de résolution ; export seul sans fenêtre interactive.
- Liens entre positions WT/mutant, connecteurs de paires de bases, filtres de
  distance et tooltips compacts ou détaillés.
- Correction du lancement des modes `minimal` et `molecule` sans métrique de
  liens explicite.
- Compatibilité des légendes et tooltips avec les versions récentes de VTK.
- Correction des colonnes PDB exportées : l'identifiant de chaîne est préservé
  à la relecture et les couleurs du viewer correspondent aux scores des rapports.

### Maintenance et compatibilité

- Tests d'intégration hors réseau et workflow de tests sous Linux et Windows.
- Métadonnées du paquet : auteur **DeeW69**, URL du dépôt corrigée et découverte
  des paquets limitée à `deepstructgenomics`.
- Les positions des substitutions dans le **Markdown** sont désormais comptées
  depuis 1, comme les hotspots. Les champs JSON existants
  `variant_analysis.substitutions[].position` et `structure.base_pairs` restent
  comptés depuis 0 pour préserver la compatibilité.
- Les scores restent des heuristiques de composition et d'appariement ; les
  projections 3D ne sont pas des structures atomiques prédites.

## 0.1.0

- Pipeline initial : séquence NCBI ou fournie, structure secondaire simplifiée,
  comparaison de variants, rapports JSON/Markdown et visualisation 3D.
