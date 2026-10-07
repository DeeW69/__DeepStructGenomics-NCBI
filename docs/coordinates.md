# Correspondance et numérotation des positions

Ce document détaille les conventions de numérotation et de repérage des positions
utilisées à travers l'ensemble de la chaîne de traitement de DeepStructGenomics-NCBI :
des séquences biologiques brutes jusqu'aux rapports d'analyse et aux fichiers de
visualisation 3D.

---

## 1. Vue d'ensemble des conventions

Deux systèmes de coordonnées coexistent dans le projet :

1. **Numérotation biologique (base 1, *1-based*)** :
   - Convention utilisée pour les positions destinées à la lecture humaine dans ce projet.
   - Le premier nucléotide d'une séquence porte le numéro **1** et le dernier le numéro **$L$** (pour une séquence de longueur $L$).
   - Utilisée pour l'ensemble des sorties destinées à la lecture humaine (rapports Markdown, messages d'erreur de validation, tooltips du viewer 3D), ainsi que dans les fichiers PDB et les fichiers de scores.

2. **Indexation interne machine (base 0, *0-based*)** :
   - Convention native des tableaux et chaînes de caractères en Python.
   - Les indices parcourent l'intervalle **$[0, L - 1]$**.
   - Présente dans les algorithmes internes (`SecondaryStructureResult`, `compare_sequences`) et conservée dans certains champs historiques du rapport JSON afin de garantir la rétrocompatibilité avec la version `v0.1.0`.

Ces positions sont **relatives à la séquence fournie**, même lorsqu'elle est
récupérée via une accession NCBI. Le pipeline ne les convertit pas en coordonnées
chromosomiques. La conversion ADN → ARN (`T` → `U`) conserve la longueur et
l'ordre des bases. Les espaces et les en-têtes FASTA ne sont pas comptés.

Dans les tableaux ci-dessous, $L$ désigne la longueur de la séquence concernée.
Le [cas des longueurs inégales](#cas-des-séquences-de-longueurs-inégales) distingue
les positions comparées dans les variants de celles utilisées pour les deltas.

---

## 2. Tableau comparatif par composant

| Composant / Fichier | Champ ou Entité | Base | Plage d'indices | Exemple concret | Rôle et justification |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Validation d'entrée**<br>`ncbi_client.py` | Message d'erreur de caractère invalide | **1** | $1 \dots L$ | `position 4` | Lisibilité biologique : permet à l'utilisateur de localiser directement l'erreur dans son fichier source (sans compter les espaces). |
| **Structure secondaire**<br>`secondary_structure.py` | Accès Python à `SecondaryStructureResult.dot_bracket[idx]` | **0** | $0 \dots L-1$ | `dot_bracket[0]` | La chaîne dot-bracket est alignée sur la séquence ; elle ne contient pas elle-même de numéros. |
| **Structure secondaire**<br>`secondary_structure.py` | `SecondaryStructureResult.base_pairs` | **0** | $0 \dots L-1$ | `(0, 11)` | Paires prédites par Nussinov sous forme de tuples d'indices Python `(i, j)` avec $i < j$. |
| **Analyse de variants**<br>`variant_analysis.py` | `VariantImpactResult.substitutions[].position` | **0** | $0 \dots L-1$ | `11` | Indice séquentiel Python de la substitution détectée. |
| **Rapport JSON**<br>`<run_id>.json` | `structure.base_pairs` | **0** | $0 \dots L-1$ | `[[0, 11]]` | **Historique v0.1.0** : paires 0-based conservées pour la compatibilité des scripts aval. |
| **Rapport JSON**<br>`<run_id>.json` | `variant_analysis.substitutions[].position` | **0** | $0 \dots L-1$ | `11` | **Historique v0.1.0** : position 0-based de la substitution. |
| **Rapport JSON**<br>`<run_id>.json` | `impact_summary.hotspots[].position` | **1** | $1 \dots L$ | `12` | **Enrichissement v0.2.0** : position 1-based identique au tableau Markdown et aux PDB. |
| **Rapport JSON**<br>`<run_id>.json` | `impact_summary.base_pair_summary` | Sans objet | Comptes et ratios | `total_pairs: 4` | Ne contient aucune position ni liste de paires ; le calcul utilise des paires converties en base 1. |
| **Rapport Markdown**<br>`<run_id>.md` | Tableau des substitutions | **1** | $1 \dots L$ | `\| 12 \| C \| A \|` | Lisibilité biologique humaine : colonne explicite `Position (1-based)`. |
| **Rapport Markdown**<br>`<run_id>.md` | Tableau des hotspots | **1** | $1 \dots L$ | `\| 12 \| -0.700 \|` | Classement des positions les plus impactées (numérotées depuis 1). |
| **Export CSV des hotspots**<br>`<run_id>_hotspots.csv`<br>`hotspots.csv` | Colonne `position` | **1** | $1 \dots L$ | `12,C,A,-0.7,0.7` | Format tabulaire pour tableurs et notebooks : position 1-based, bases WT/mutant et deltas. |
| **PDB exportés par le pipeline**<br>`wt_structure.pdb`<br>`mutant_structure.pdb` | Colonne numéro de résidu `resSeq` | **1** | $1 \dots L$ | Résidu `A:12:` | Un pseudo-atome par nucléotide ; numérotation continue depuis 1 sur la chaîne `A`. |
| **Fichiers de scores**<br>`wt_scores.json`<br>`mutant_scores.json` | `position_scores` (clés JSON) | **1** | `"1"` à `str(L)` | `"12": 0.3` | Clés textuelles correspondant au numéro de résidu `resSeq` du PDB exporté. |
| **Fichiers de scores**<br>`wt_scores.json`<br>`mutant_scores.json` | `residue_scores` (clés compactes) | **1** | Résidus 1 à $L$, chaîne `A` | `"A:12:": 0.3` | Format du projet `chaîne:resSeq:icode`, produit par `ResidueKey.as_compact()`. |
| **Résumé d'impact**<br>`impact_summary.json` | `hotspots[].position` | **1** | $1 \dots L$ | `12` | Fichier autonome copié dans le dossier `visualization/`, identique à la clé du rapport JSON. |
| **Connecteurs 3D**<br>`--base-pairs-wt`<br>`--base-pairs-mut` | Liste de paires JSON `base_pairs` | **1** | $1 \dots L$ | `[[1, 12]]` | Paires 1-based fournies au visualiseur VTK ; les paires contenant un indice $\le 0$ sont ignorées. |
| **Viewer interactif**<br>`scripts/view_3d.py` | Tooltips des overlays | **1** | $1 \dots L$ pour les exports du pipeline | `pos=12`, `residue=A:12:` | La position séquentielle et la clé de résidu sont affichées séparément. |

---

## 3. Détail étape par étape

Les nouveaux exports sur `main` ajoutent `secondary_structures.json` : les
objets `reference` et `mutant` contiennent `sequence`, `dot_bracket`, `base_pairs`
(indices depuis **0**) et `scores` (liste alignée sur la séquence, accès Python
depuis **0**). `mutant` vaut `null` sans comparaison. Le viewer secondaire ajoute
1 aux indices pour l'affichage. Les fichiers séparés `wt_base_pairs.json` et
`mut_base_pairs.json` contiennent déjà les paires depuis **1** pour VTK ; ne pas
leur ajouter 1 une seconde fois.

### 3.1. Entrées et validation des séquences

Lors de l'ingestion d'une séquence (accès NCBI, fichier FASTA `--fasta`, ou chaîne CLI `--sequence`) :

- Les espaces et sauts de ligne sont d'abord éliminés : la séquence compacte résultante a une longueur $L$.
- La fonction `normalize_rna_sequence` vérifie chaque caractère. Si un caractère interdit est rencontré (code dégénéré, gap, chiffre, symbole) :
  ```text
  Caractere non pris en charge 'N' a la position 4 (numerotation depuis 1, sans espaces). Bases acceptees : A, C, G, T, U.
  ```
- La position rapportée commence à **1**. L'utilisateur peut ainsi inspecter directement son fichier FASTA au 4ᵉ nucléotide.

### 3.2. Traitement interne (Python)

Les structures internes utilisent l'indexation Python native :

- **Dot-bracket** : la chaîne commence à l'index `0`. Le premier nucléotide est `seq[0]`, son état de repliement est `dot_bracket[0]`.
- **Appariements** : `base_pairs` est une liste de tuples `(i, j)` où $0 \le i < j \le L - 1$.
- **Substitutions** : la fonction `compare_sequences` compare `reference_seq[idx]` et `mutant_seq[idx]`. Le champ `position` dans chaque entrée de substitution vaut `idx` (entre $0$ et $\max(L_{\text{wt}}, L_{\text{mut}}) - 1$).

### 3.3. Génération des rapports (Markdown et JSON)

Lors de la sérialisation via `report_generator.py` :

- **Dans le fichier Markdown (`<run_id>.md`)** :
  - La table des substitutions convertit l'indice en ajoutant `+ 1` :
    ```markdown
    | Position (1-based) | Reference | Mutant |
    | 12 | C | A |
    ```
  - La section des variations par position précise explicitement :
    > « Positions numérotées à partir de 1. Delta = score mutant - score référence. »
  - Les hotspots affichent directement cette position 1-based.
- **Dans le fichier JSON (`<run_id>.json`)** :
  - `structure.base_pairs` contient les paires 0-based d'origine (ex. `[[0, 11]]`).
  - `variant_analysis.substitutions[].position` contient l'index 0-based d'origine (ex. `11`).
  - `impact_summary.hotspots[].position` contient la position 1-based (ex. `12`).
  - Cette distinction préserve les automatisations écrites pour la version 0.1.0 tout en fournissant des métadonnées 1-based alignées avec le reste des artefacts.

### 3.4. Artefacts de visualisation 3D

Le dossier `outputs/<run_id>/visualization/` contient les fichiers consommés par `scripts/view_3d.py` :

- **Fichiers PDB (`wt_structure.pdb`, `mutant_structure.pdb`)** :
  - Dans les exports du pipeline, chaque nucléotide est modélisé par un pseudo-atome `P` du résidu `NTP` sur la chaîne `A`.
  - La numérotation atomique `serial` et la numérotation résiduelle `resSeq` commencent toutes deux à **1** :
    ```text
    ATOM      1  P   NTP A   1      12.000   0.000   0.000  1.00  0.00           P
    ...
    ATOM     12  P   NTP A  12      -2.341 -11.769  30.800  1.00  0.00           P
    ```
- **Fichiers de scores (`wt_scores.json`, `mutant_scores.json`)** :
  - `derive_position_scores_from_structure` convertit l'index `idx` en `idx + 1` :
    les clés suivantes illustrent un extrait de `mutant_scores.json` de la démo.
    ```json
    {
      "version": 1,
      "residue_scores": {
        "A:1:": 0.40,
        "A:12:": 0.30
      },
      "position_scores": {
        "1": 0.40,
        "12": 0.30
      }
    }
    ```
  - Les clés de `residue_scores` (`"A:1:"`, `"A:12:"`) sont construites par `ResidueKey.as_compact()` au format `<chaîne>:<resSeq>:<icode>` ; le champ final est vide en l'absence de code d'insertion. Ce format appartient au projet. Le lecteur utilise Biopython pour extraire ces éléments du PDB.
  - Après `json.load`, les clés de `position_scores` sont des chaînes : `payload["position_scores"]["12"]`. Après `load_score_table`, elles sont converties en entiers : `table.position_scores[12]`. Ces accès ne sont pas interchangeables.
- **Connecteurs de paires de bases (`--base-pairs-wt`, `--base-pairs-mut`)** :
  - Doivent contenir des paires **1-based** : `{"base_pairs": [[1, 12]]}`.
  - La fonction `build_base_pair_segments` soustrait 1 (`coords[i - 1]`) pour retrouver les coordonnées NumPy correspondantes.
  - Le chargeur ignore les paires avec un indice nul ou négatif, sans les convertir. Copier directement les paires 0-based du rapport perdrait donc certaines paires et décalerait les autres. Ajouter `1` aux **deux** extrémités avant de créer le fichier de connecteurs.

**Structures importées :** un PDB ou mmCIF externe conserve ses identifiants de
chaîne et de résidu ; le lecteur ne les renumérote pas. Dans les tooltips des
overlays, la position utilise `sequence_index` si disponible, sinon le rang
de l'atome lu, compté depuis 1. L'égalité « position = `resSeq` = numéro d'atome »
concerne les exports du pipeline à un pseudo-atome par nucléotide ; elle ne
doit pas être supposée pour une structure externe avec plusieurs atomes par résidu.

---

## 4. Exemple fil conducteur (démo synthétique 12 nt)

Soit les deux ARN synthétiques de la démo hors ligne (`data/examples/rna_variant.json`) :

- **WT** : `G G G G A A A A C C C C`
- **Mutant** : `G G G G A A A A C C C A` (substitution C $\rightarrow$ A en dernière position)

### Tableau de correspondance position par position

| Base WT | Base MUT | Index Python (0-based) | Position biologique (1-based) | Clé PDB (`residue_key`) | Score WT | Score MUT | Delta (MUT - WT) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **G** | **G** | `0` | **1** | `A:1:` | 1.000 | 0.400 | **-0.600** |
| **G** | **G** | `1` | **2** | `A:2:` | 1.000 | 1.000 | 0.000 |
| **G** | **G** | `2` | **3** | `A:3:` | 1.000 | 1.000 | 0.000 |
| **G** | **G** | `3` | **4** | `A:4:` | 1.000 | 1.000 | 0.000 |
| **A** | **A** | `4` | **5** | `A:5:` | 0.300 | 0.300 | 0.000 |
| **A** | **A** | `5` | **6** | `A:6:` | 0.300 | 0.300 | 0.000 |
| **A** | **A** | `6` | **7** | `A:7:` | 0.300 | 0.300 | 0.000 |
| **A** | **A** | `7` | **8** | `A:8:` | 0.300 | 0.300 | 0.000 |
| **C** | **C** | `8` | **9** | `A:9:` | 1.000 | 1.000 | 0.000 |
| **C** | **C** | `9` | **10** | `A:10:` | 1.000 | 1.000 | 0.000 |
| **C** | **C** | `10` | **11** | `A:11:` | 1.000 | 1.000 | 0.000 |
| **C** | **A** | `11` | **12** | `A:12:` | 1.000 | 0.300 | **-0.700** |

Les scores incluent le bonus lié à la distance entre partenaires, puis sont
plafonnés à 1. Le `G` non apparié en position 1 du mutant vaut 0,4 : socle 0,3
et bonus de composition 0,1. Ces scores sont heuristiques.

### Représentation dans les artefacts produits

- **Paire terminale rompue** :
  - En interne Python et dans `payload["structure"]["base_pairs"]` : `[0, 11]`
  - Dans un fichier de connecteurs 3D préparé depuis les paires WT : `[1, 12]`.
  - Le PDB contient les résidus `A:1:` et `A:12:`, mais **pas l'appariement** entre eux. Les connecteurs se lisent dans un fichier JSON séparé.
- **Substitution C $\rightarrow$ A** :
  - Dans `payload["variant_analysis"]["substitutions"]` : `{"position": 11, "reference": "C", "mutant": "A"}`
  - Dans le tableau Markdown : `| 12 | C | A |`
- **Hotspots** :
  - Dans `impact_summary.json` et le rapport Markdown : position `12` ($\Delta = -0.700$) puis position `1` ($\Delta = -0.600$).

---

## 5. Formules de conversion et recommandations pour les scripts aval

Pour consommer les résultats du pipeline sans ambiguïté :

Lancer d'abord `python scripts/run_demo.py` depuis la racine du dépôt. Les
exemples ci-dessous lisent le rapport produit dans `outputs_demo/`.

### En Python

```python
import json
from pathlib import Path

with open("outputs_demo/offline_demo.json", encoding="utf-8") as f:
    report = json.load(f)

# Accès 0-based (historique v0.1.0)
for sub in report["variant_analysis"]["substitutions"]:
    idx_0based = sub["position"]
    pos_1based = idx_0based + 1
    print(f"Index Python: {idx_0based} -> Position biologique: {pos_1based}")

# Accès 1-based direct (enrichissement v0.2.0)
for hotspot in report["impact_summary"]["hotspots"]:
    pos_1based = hotspot["position"]
    idx_0based = pos_1based - 1
    print(f"Hotspot à la position {pos_1based} (delta = {hotspot['delta']})")

# Conversion explicite des paires WT pour --base-pairs-wt.
wt_pairs = [[i + 1, j + 1] for i, j in report["structure"]["base_pairs"]]
Path("outputs_demo/wt_base_pairs.json").write_text(
    json.dumps({"base_pairs": wt_pairs}, indent=2), encoding="utf-8"
)
```

### En R

Cet exemple nécessite le paquet R `jsonlite`, distinct des dépendances Python.

```R
library(jsonlite)
report <- fromJSON("outputs_demo/offline_demo.json")

# R indexe nativement en base 1
hotspots_pos <- report$impact_summary$hotspots$position # direct: 12, 1
substitutions_pos <- report$variant_analysis$substitutions$position + 1 # conversion en base 1
```

### Cas des séquences de longueurs inégales

La comparaison est positionnelle, **sans alignement**. Le tableau des
substitutions parcourt les indices $0 \dots \max(L_{\text{wt}}, L_{\text{mut}})-1$ ;
au-delà de la séquence la plus courte, `reference` ou `mutant` vaut `"-"` pour
signaler une base absente. Ce remplissage ne constitue pas un alignement et
ne localise pas une insertion ou une délétion biologique.

Les deltas de scores et les hotspots portent uniquement sur les positions
partagées $1 \dots \min(L_{\text{wt}}, L_{\text{mut}})$. Ainsi, avec `AUGC` et
`AUGCA`, le rapport des variants inclut `{"position": 4, "reference": "-", "mutant": "A"}`
(position 5 dans le Markdown), mais aucun delta de score n'est calculé à cette position.
