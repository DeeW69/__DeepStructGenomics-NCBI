# Release Checklist

Guide concis pour préparer, taguer et publier une version de DeepStructGenomics-NCBI.

## 1. Préparer la release

1. Confirmer que l’arbre est propre :
   ```powershell
   git status -sb
   ```
2. Installer le verrou dans un environnement neuf selon le
   [guide des dépendances](dependencies.md), puis vérifier la compatibilité et
   lancer les tests avec le Python de cet environnement :
   ```powershell
   python -m pip check
   python -m pytest
   ```
3. Vérifier que `scripts/view_3d.py` et `docs/visualization_3d.md` sont alignés (options/behaviour).
   Relire rapidement les exemples, surtout si des flags ont évolué.
4. Générer un run de démonstration complet pour valider les artefacts :
   ```powershell
   python scripts/run_demo.py
   python scripts/plot_demo.py
   ```
   - Inspecter `outputs_demo/offline_demo.md` et `outputs_demo/offline_demo/visualization/`
     (scores, manifest, impact_summary), puis comparer avec `docs/demo.md`.
   - Capturer au moins un PNG (ex. `delta.png`) via `scripts/view_3d.py --export`.
   - Vérifier que la plage et les couleurs du viewer correspondent aux deltas
     du rapport (pour la démo : positions 1 = -0,6 et 12 = -0,7).
5. Vérifier le paquet installable, si l'outil `build` est installé :
   ```powershell
   python -m build
   ```

## 2. Versioning

- Convention : `vMAJOR.MINOR.PATCH` (ex. `v0.1.0`, `v0.2.0`, `v1.0.0`).
- Incréments recommandés :
  - **PATCH** : corrections mineures, docs, refactors sans impact API.
  - **MINOR** : nouvelles fonctionnalités rétro-compatibles, nouveaux artefacts.
  - **MAJOR** : changements incompatibles, restructuration CLI/API, réécriture pipeline.

## 3. Créer un tag

1. Tag local annoté :
   ```bash
   git tag -a vX.Y.Z -m "Release vX.Y.Z"
   ```
2. Vérifier les tags :
   ```bash
   git tag --list
   ```
3. Pousser le tag :
   ```bash
   git push origin vX.Y.Z
   ```

## 4. Release GitHub

Vérifier que les deux jobs Linux/Windows réussissent pour le commit à publier.

1. Depuis l’onglet **Releases**, cliquer « Draft a new release ».
2. Sélectionner le tag `vX.Y.Z`, donner un titre clair (`DeepStructGenomics vX.Y.Z`).
3. Rédiger les notes :
   - Highlights / nouvelles features.
   - Bug fixes.
   - Breaking changes / migrations (si applicable).
4. Joindre les artefacts utiles (ex. `overlay.png`, `delta.png`, un `visualization_manifest.json` d’exemple).
5. Publier (ou marquer en pré-release si alpha/beta, cf. section suivante).

## 5. Pre-release checklist (alpha/beta)

- Mention explicite dans le tag (`v0.x.y-alpha` ou `v1.0.0-beta.1`).
- Notifier que les artefacts peuvent changer et que l’API/CLI n’est pas figée.
- Documenter les éléments à valider avant GA :
  - Couverture tests minimale.
  - Feedback sur CLI / docs.
  - Revue manuelle des artefacts (PNG, impact_summary, manifest).
- Utiliser l’option « This is a pre-release » lors de la publication GitHub.

## 6. Rollback / suppression de tag

> ⚠️ À n’utiliser que si la release contient un problème bloquant.

1. Supprimer le tag local :
   ```bash
   git tag -d vX.Y.Z
   ```
2. Supprimer le tag distant :
   ```bash
   git push origin :refs/tags/vX.Y.Z
   ```
3. Si la release GitHub existe, la supprimer également (UI GitHub) et communiquer clairement sur le rollback.
