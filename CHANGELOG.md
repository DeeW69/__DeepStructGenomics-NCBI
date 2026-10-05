# Historique des versions

## 0.2.0 — en préparation

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
