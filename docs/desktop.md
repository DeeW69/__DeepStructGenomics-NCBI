# Espace d'analyse desktop

Cette première interface prépare **v0.5.0** sur `main`. La dernière release
reste v0.4.6. Elle utilise PySide6/Qt et appelle les fonctions scientifiques
existantes ; les viewers et commandes CLI restent disponibles.

## Installation et lancement

Depuis un checkout avec Python 3.10 ou ultérieur :

```bash
python -m pip install --require-hashes --only-binary=:all: -r requirements-lock.txt -r requirements-gui-lock.txt
# Facultatif pour ARN, nécessaire pour les calculs Hi-C avancés
python -m pip install --require-hashes --only-binary=:all: -r requirements-research-lock.txt
python scripts/run_gui.py
```

Sous Windows, remplacer `python` par le chemin de l'interpréteur de votre
environnement, par exemple `.\.venv\Scripts\python.exe`.

```bash
python scripts/run_gui.py --workspace outputs_mes_analyses
python scripts/run_gui.py --manifest outputs_demo/offline_demo/visualization/visualization_manifest.json
```

Après installation du wheel avec l'extra `gui`, utiliser `deepstructgenomics`
ou `python -m deepstructgenomics.gui.app` depuis n'importe quel répertoire.
Les fichiers d'exemple Hi-C sont fournis dans le dépôt et l'archive source.

## Parcours ARN

1. Cliquer sur **Essayer la démo C12A**, ou ouvrir **Séquences ARN**.
2. Choisir une saisie directe, un FASTA mono-entrée ou une accession NCBI
   nuccore exacte. La recherche par nom de gène n'est pas encore disponible.
3. Ajouter éventuellement une séquence ou un FASTA mutant, puis lancer l'analyse.
4. La page **Comparaison** affiche les paires et les différences de scores.
   Cliquer sur une position du graphique, une ligne du tableau ou utiliser
   le sélecteur pour voir bases, scores et partenaires d'appariement.
5. L'onglet **Vue 3D illustrative** charge la scène VTK existante à la demande.
   Référence, mutant, squelette et paires peuvent être masqués séparément.
   Un clic sur un résidu mutant renseigne l'inspecteur. Les paramètres internes
   sont dans **Informations techniques**.

La démo synthétique contient 12 bases : WT 4 paires, mutant 3, paire 1–12 perdue,
C12A, Δ maximal absolu 0,700. Le graphique reprend les résultats exportés,
sans nouvelle prédiction. Nussinov pondéré ne fournit pas d'énergie MFE.
Les séquences de longueurs différentes sont comparées par position, sans alignement.
Une valeur absente apparaît comme indisponible, jamais comme un zéro mesuré.

La géométrie 3D ARN dépend de la longueur. Elle ne montre pas un repliement
moléculaire ni un déplacement physique causé par la mutation. La vue superposée
requiert un mutant ; une référence seule reste consultable en 2D et dans le tableau.

## Parcours Hi-C

Choisir un TSV de comptes bruts, renseigner assemblage, chromosome, début,
fin exclue et taille de bin. Déclarer explicitement que les contacts absents
sont des zéros observés avant de lancer le calcul. **Remplir l'exemple synthétique**
charge les paramètres des données locales, sans lancer automatiquement l'analyse.

La région est limitée à 400 bins ; les coordonnées sont en bp, base 0.
Cette première interface conserve les paramètres par défaut du moteur : fenêtre
3 bins, 999 permutations, graine 46, FDR 0,05. Les paramètres supplémentaires
restent accessibles par CLI. Voir [les méthodes et limites](hic_analysis.md).

La liste de vues propose les contacts équilibrés, les boucles candidates et
la reconstruction MDS. Les lignes cyan indiquent les frontières retenues.
Une non-convergence est signalée et ne produit ni tests ni reconstruction.
Les q-values dépendent des modèles nuls ; les coordonnées sont inférées en unités
arbitraires. Le caractère synthétique de l'exemple ne constitue pas une validation biologique.

## Résultats et calculs

Les calculs s'exécutent dans un processus Python séparé. L'interface reste
disponible pour parcourir les résultats précédents ; un seul calcul tourne à la
fois. **Annuler le calcul** arrête le processus. Les éventuels fichiers partiels
restent dans son dossier propre, sans entrée ajoutée à l'historique. Attendre
la fin du calcul ou l'annuler avant de fermer l'application.

Chaque lancement crée un dossier unique dans `outputs_gui` (ou `--workspace`).
Les exports sont ceux du moteur : JSON/Markdown/CSV pour ARN, JSON/Markdown/TSV/BED
pour Hi-C. Les boutons PNG exportent la figure affichée. **Dossier des résultats**
ouvre les fichiers disponibles ; **Ouvrir un résultat** accepte un manifeste ARN
ou un `hic_analysis.json`. Les anciens manifests dépourvus de structures secondaires
doivent être régénérés par le pipeline actuel.

L'historique local `recent.json` conserve les chemins des 20 derniers résultats,
pas les séquences saisies ni l'e-mail NCBI. Un résultat déplacé doit être rouvert
à son nouvel emplacement. Les fichiers exportés et le cache contiennent les données
analysées. Aucun transfert réseau n'a lieu en dehors des demandes NCBI.

## Architecture et validation

- `gui/services.py` délègue aux pipelines, lit les artefacts et gère l'historique.
- `gui/worker.py` reçoit la demande sur stdin et renvoie le chemin des résultats.
- `gui/main_window.py` supervise le processus et les quatre pages.
- `gui/comparison.py`, `gui/hic.py` et `gui/vtk_viewer.py` affichent les résultats.

L'intégration VTK utilise le widget officiel
[QVTKRenderWindowInteractor](https://docs.vtk.org/en/v9.6.1/api/python/vtkmodules/vtkmodules.qt.QVTKRenderWindowInteractor.html).
Le thème est empaqueté avec l'application ; aucune interface web ni service
d'hébergement n'est nécessaire.

Les tests couvrent la démo, les erreurs suivies d'une nouvelle analyse, l'annulation,
la réouverture, les mutants absents/de longueurs différentes et le parcours Hi-C.
Le rendu OpenGL dépend du poste ; en cas d'indisponibilité, la comparaison 2D
reste accessible. Recherche NCBI avec aperçu, lots FASTA, prédicteur ViennaRNA
dans l'interface et gestion de projets multiples sont des étapes ultérieures.
