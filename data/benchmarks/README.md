# Références ARN expérimentales

`bprna_pdb.json` contient 30 séquences et annotations secondaires issues du
sous-ensemble **PDB** de bpRNA-1m v1.0. Les références sont dérivées de structures
expérimentales déposées ; elles ne sont pas produites par nos prédicteurs.

Source : [bpRNA-1m](https://bprna.cgrb.oregonstate.edu/about.php),
[article de Danaee et al.](https://doi.org/10.1093/nar/gky285),
[archive ST](https://bprna.cgrb.oregonstate.edu/bpRNA_1m/stFiles.zip).
Ces données conservent leur attribution à leurs producteurs ; la licence MIT
du logiciel ne constitue pas une nouvelle licence des données sources.

La sélection précède toute prédiction : 669 annotations PDB examinées, longueur
20–120 nt, alphabet ACGU, notation `.()` sans pseudonœuds, séquences identiques
dédupliquées. Parmi 207 candidats, tirage de 30 entrées avec graine 45.
Les SHA-256 de l'archive complète et de chaque membre sont enregistrés.
L'archive originale volumineuse n'est pas incluse dans le dépôt.

Après téléchargement de l'archive, reproduire le corpus :

```bash
python scripts/prepare_rna_benchmark.py --archive data/raw/bprna/stFiles.zip --output outputs_benchmark/references.json
```

Limites : petit échantillon, absence de filtrage des homologues et d'équilibrage
par famille, structures sans pseudonœuds uniquement. Les contacts non canoniques
et boucles courtes des références sont conservés, même si les modèles ne peuvent
pas tous les reproduire. Le contexte cristallin, les protéines et ligands ne sont
pas modélisés. L'indépendance des références vis-à-vis du projet ne signifie pas
qu'elles sont absentes de l'ajustement historique des paramètres de ViennaRNA.
