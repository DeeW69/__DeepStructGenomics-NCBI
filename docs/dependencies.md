# Dépendances et environnement reproductible

`requirements-lock.txt` fixe les versions des dépendances du pipeline et des
tests, y compris leurs dépendances transitives. Les empreintes SHA-256 permettent
à `pip` de vérifier les distributions téléchargées. Ce fichier est généré ; ne
pas modifier ses versions ou ses empreintes à la main.

`pyproject.toml` décrit les dépendances du paquet avec leurs versions minimales.
`requirements.txt` contient ces mêmes dépendances ainsi que `pytest` pour le
développement. Ces fichiers conservent des plages de versions ; le verrou est
l'environnement de référence utilisé par la CI et les commandes ci-dessous.

## Installer l'environnement de référence

Créer un environnement virtuel neuf avec Python 3.10 ou ultérieur.

Linux / macOS :

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --require-hashes --only-binary=:all: -r requirements-lock.txt
python -m pip check
python -m pytest
```

Windows / PowerShell, sans activation :

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install --require-hashes --only-binary=:all: -r requirements-lock.txt
.\.venv\Scripts\python.exe -m pip check
.\.venv\Scripts\python.exe -m pytest
```

Les commandes utilisent des roues précompilées (`--only-binary=:all:`) afin
d'éviter une compilation locale avec des dépendances de construction non
verrouillées. La CI valide Linux avec Python 3.10 et Windows avec Python 3.12.
La résolution universelle contient les marqueurs Python et système nécessaires
à ces environnements ; elle ne garantit pas la disponibilité de roues pour
toutes les architectures ou les futures versions de Python.

Le verrou couvre l'exécution et les tests. Il ne verrouille pas l'interpréteur,
les pilotes graphiques, les bibliothèques système, ni les outils de construction
du paquet (`build`, `setuptools`, `wheel`). Une installation ordinaire du paquet
avec `pip install .` utilise les plages de `pyproject.toml` et ne consomme pas
automatiquement ce verrou.

## Régénérer le verrou

### Interface desktop (préparation v0.5.0)

PySide6 est facultatif ; le socle CLI ne l'importe pas. Après le verrou principal :

```bash
python -m pip install --require-hashes --only-binary=:all: -r requirements-gui-lock.txt
python scripts/run_gui.py
```

Le paquet expose l'extra `gui` et le lanceur `deepstructgenomics`. Le thème Qt est
inclus dans le wheel. Les méthodes Hi-C avancées nécessitent aussi `research`.
Le verrou GUI fixe PySide6 et ses dépendances Qt/Shiboken ; il se régénère avec :

```bash
uv pip compile requirements-gui.txt --constraint requirements-lock.txt --universal --python-version 3.10 --generate-hashes --output-file requirements-gui-lock.txt --no-python-downloads
```

La CI exécute les parcours desktop avec `QT_QPA_PLATFORM=offscreen` après les
tests du socle et de la recherche. Le rendu VTK/OpenGL interactif doit aussi être
vérifié sur le poste cible ; les tests hors écran ne valident pas les pilotes.

### Extensions de recherche (depuis v0.4.5)

Après le verrou principal, installer les prédicteurs et statistiques facultatifs :

```bash
python -m pip install --require-hashes --only-binary=:all: -r requirements-research-lock.txt
```

ViennaRNA 2.7.2 exécute la prédiction MFE. SciPy fournit les probabilités statistiques
de l'analyse Hi-C avancée : 1.15.3 sous Python 3.10, 1.17.1 à partir de Python 3.11
pour rester compatible avec les versions NumPy du verrou principal.
Ces versions sont aussi déclarées dans l'extra `research` du paquet. La CI teste
le socle avant leur installation, puis l'ensemble avec ces extensions.

Le verrou optionnel se régénère en conservant les contraintes du socle :

```bash
uv pip compile requirements-research.txt --constraint requirements-lock.txt --universal --python-version 3.10 --generate-hashes --output-file requirements-research-lock.txt --no-python-downloads
```

### Verrou principal

La version du générateur utilisée est **uv 0.12.5**. `uv` est nécessaire seulement
pour régénérer le fichier, pas pour installer le projet ou lancer les tests.
Si besoin, l'installer dans un environnement de maintenance distinct :

```bash
python -m pip install uv==0.12.5
```

Depuis la racine du dépôt :

```bash
uv pip compile requirements.txt --universal --python-version 3.10 --generate-hashes --output-file requirements-lock.txt
```

`--universal` conserve les marqueurs de plateforme et de version Python ;
`--python-version 3.10` impose ici la borne minimale Python 3.10, même si la
commande est lancée avec un Python plus récent. Sans option de mise à jour,
`uv` privilégie les versions déjà présentes dans le verrou lorsque celles-ci
respectent toujours les contraintes.

Pour mettre à jour une dépendance, ajouter `--upgrade-package NOM` à cette
commande ; pour réévaluer toutes les versions, ajouter `--upgrade`. Lorsqu'une
dépendance du pipeline est ajoutée ou modifiée, aligner `pyproject.toml` et
`requirements.txt` avant la régénération.

Relire le diff du verrou, réinstaller dans un environnement neuf avec la commande
à empreintes ci-dessus, puis lancer `pip check`, les tests et la démo hors ligne
(`python scripts/run_demo.py`). Attendre les deux jobs Linux/Windows avant de
publier une nouvelle release fondée sur ce verrou. Les changements du fichier
source et du verrou doivent être livrés dans le même commit.
