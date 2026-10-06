# Historique des versions

## Non publié — stabilisation v0.2.x

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
