# Configuration de l'alignement — 0.5.0rc1

Une seule classe indépendante de Qt, `AlignmentConfig`, alimente le moteur,
le pipeline, les commandes et le worker GUI. L'API Python pure utilise ses
valeurs par défaut ; les points d'entrée chargent la configuration externe puis
transmettent un instantané à l'analyse. Aucun réglage de `pairwise.py` n'est à éditer.

## Priorité et fichier externe

**Défauts < `.env` < environnement système < arguments CLI explicites < choix GUI
pour l'analyse courante.** Un paramètre CLI omis ne masque pas l'environnement.
La GUI commence avec la configuration CLI/environnement et ne modifie ni `.env`
ni les paramètres des analyses antérieures.

Le fichier automatiquement recherché est `.env` dans le **répertoire courant**.
`--env-file chemin/config.env` sélectionne explicitement un fichier externe ; un
chemin explicite introuvable est une erreur. Sans `.env`, les défauts suffisent.
Le parseur léger n'ajoute aucune dépendance : lignes `DSG_NOM=valeur`, lignes
vides et commentaires sur une ligne commençant par `#`, guillemets simples/doubles
facultatifs autour de toute la valeur. Pas d'expansion de variables, pas de commande,
pas de commentaire en fin de valeur. Les autres préfixes sont ignorés. Les chemins
FASTA sont relatifs au répertoire courant, pas au fichier `.env`.

Les limites doivent être des entiers strictement positifs, les scores des nombres
finis, `enabled` vaut `true/false` ou `1/0`, et la politique vaut `reject/positional`.
Une valeur invalide produit une erreur explicite, sans retour silencieux au défaut.

Copier au besoin [`.env.example`](../.env.example) vers `.env` : ce dernier reste
ignoré par Git. Aucun des deux fichiers n'est embarqué dans le wheel ou le sdist.

```dotenv
DSG_ALIGNMENT_ENABLED=true
DSG_MAX_ALIGNMENT_LENGTH=10000
DSG_MAX_ALIGNMENT_CELLS=4000000
DSG_ALIGNMENT_MATCH=2
DSG_ALIGNMENT_MISMATCH=-1
DSG_ALIGNMENT_GAP_OPEN=-5
DSG_ALIGNMENT_GAP_EXTEND=-1
DSG_ALIGNMENT_OVERFLOW_POLICY=reject
# Une source par côté : séquence OU FASTA
# DSG_WT_SEQUENCE=ACGU
# DSG_MUT_SEQUENCE=ACGGU
# DSG_WT_FASTA=data/wt.fasta
# DSG_MUT_FASTA=data/mut.fasta
```

Une source explicite remplace la source du même côté venant d'une couche inférieure.
Fournir deux sources WT (ou MUT) dans une même couche est une erreur. Un FASTA doit
contenir une seule entrée ; la CLI de pipeline conserve son mode lot séparé.

## Commandes PowerShell

Après installation des dépendances, depuis le dépôt :

```powershell
python scripts/run_alignment.py --show-config
python scripts/run_alignment.py --wt ACGU --mut ACGGU
python scripts/run_alignment.py --wt-fasta wt.fasta --mut-fasta mut.fasta
python scripts/run_alignment.py --env-file C:/mes-reglages/dsg.env --max-cells 25000000 --show-config
python scripts/run_pipeline.py --sequence GGGGAAAACCCC --mutant-sequence GGGGUAAAACCCA --max-cells 4000000 --output-dir outputs_indel
python scripts/run_gui.py --max-cells 25000000 --sequence ACGU --mutant-sequence ACGGU
```

Les options communes sont `--alignment/--no-alignment`, `--max-length`,
`--max-cells`, `--match`, `--mismatch`, `--gap-open`, `--gap-extend`,
`--overflow-policy`. Les noms historiques `--sequence`, `--mutant-sequence`,
`--fasta`, `--mutant-fasta` restent disponibles dans le pipeline et la GUI.
Après installation du wheel, `dsg-align` et `dsg-benchmark-alignment` remplacent
respectivement `python scripts/run_alignment.py` et `python scripts/benchmark_alignment.py`,
depuis n'importe quel dossier.

Une variable PowerShell reste possible, sans être obligatoire :

```powershell
$env:DSG_MAX_ALIGNMENT_CELLS="50000000"
python scripts/run_alignment.py --show-config
Remove-Item Env:DSG_MAX_ALIGNMENT_CELLS
```

## GUI : paramètres avancés

Dans **Séquences → Paramètres avancés**, l'alignement peut être activé/désactivé.
Les profils Standard **4 M**, Étendu **25 M**, Très étendu **50 M** changent la limite
et rétablissent le scoring et la longueur maximale par défaut. Le profil
**Personnalisé** rend ces champs éditables. Le bouton de restauration rétablit
tous les défauts, dont le rejet en cas de dépassement.

La limite signifie « travail maximal autorisé », pas « travail effectué » :
50 M autorise environ 7 071 × 7 071 nt, mais WT 12 × MUT 12 représente seulement
144 cellules, soit 0,000288 % de cette limite. `alignment_cell_count` définit
partout **len(WT) × len(MUT)**, sans ligne/colonne de bord ; ce n'est ni une mesure
des allocations internes de Biopython ni une estimation exacte de la RAM.

Le formulaire affiche le travail théorique des séquences directement saisies.
Pour les FASTA, les longueurs sont disponibles après lecture par le worker et
dans **Comparaison → Paramètres utilisés**, avec la configuration **enregistrée**.
Le tableau WT/MUT est virtualisé : Qt demande les cellules visibles, sans créer
tous les objets de 7 000 lignes à l'ouverture.

Au dépassement, aucun alignement n'est tenté. Les actions permettent de choisir
la comparaison positionnelle, modifier les limites ou annuler ; il faut relancer
explicitement l'analyse. Aucune limite n'est relevée automatiquement. Le mode
positionnel est signalé et peut décaler les correspondances lors d'indels.

## Benchmark synthétique, séparé du pipeline ARN

```powershell
python scripts/benchmark_alignment.py --length 7070 --max-cells 50000000
python scripts/benchmark_alignment.py --length 7070 --max-cells 50000000 --pattern deterministic-random --seed 42 --substitutions 10 --insertions 3 --deletions 3 --json-output outputs/benchmarks/50m.json
python scripts/benchmark_alignment.py --length 5000 --max-cells 25000000 --pattern repetitive
python scripts/run_alignment.py --generate-length 7070 --mutations 10 --insertions 2 --deletions 2 --seed 42 --max-cells 50000000
```

La génération est déterministe pour une graine donnée, en alphabet ACGU. Les
substitutions et délétions ciblent des positions WT distinctes ; les insertions
sont ensuite ajoutées au mutant. Les nombres de différences **retrouvés par
l'alignement** peuvent différer des opérations générées, surtout dans les répétitions.
Les données sont signalées comme synthétiques et ne passent jamais par le
repliement, les rapports biologiques ou l'historique des analyses GUI.

Le résumé affiche statut, temps, pic mémoire du processus, identité, score,
substitutions/insertions/délétions et ambiguïté ; `--verbose` seul affiche les
chaînes complètes. Un rejet retourne le code 2. La mémoire est le pic du working
set Windows / RSS Unix du **processus entier**, imports compris, en MiB : ce n'est
pas la mémoire additionnelle du moteur. Si indisponible, la mesure vaut `null`.
Lancer chaque mesure dans un processus neuf. Les grandes tailles restent hors CI.

Voir les [mesures Windows réelles](benchmarks/alignment_windows.md). Elles ne
justifient pas de relever automatiquement le défaut de production. Les limites
d'alignement ne protègent pas du coût du repliement Nussinov ; accepter 50 M de
cellules dans ce benchmark ne garantit pas une analyse ARN complète rapide.

## Traçabilité et ambiguïtés

Les blocs JSON `alignment` et `alignment_run` conservent scoring, limites,
longueurs, nombre de cellules, statut et règle de départage. `ambiguous=true`
signale plusieurs optima ; la correspondance déterministe sert à l'affichage,
sans certitude biologique sur le placement du gap. Les rapports Markdown gardent
une section technique concise. Voir [coordonnées et migration](alignment.md).
