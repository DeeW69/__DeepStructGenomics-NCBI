# Ameliorations du repliement ARN

## Perimetre
- [ ] Benchmark complet de `deepstructgenomics/rna/secondary_structure.py`
- [ ] Nouvelle energie / metrique dans la matrice de pairing
- [ ] Ajustement des sorties JSON/Markdown (`deepstructgenomics/reporting/report_generator.py`)
- [ ] Parametrage pipeline/CLI (`deepstructgenomics/pipeline.py`, `scripts/run_pipeline.py`)
- [ ] Mise a jour documentation / exemples de `data/examples/`
Indiquez quels volets sont concernes, et precisez si vous ciblez uniquement certains types d'ARN (miRNA, lncRNA, motifs viraux, etc.).

## Objective
Formuler clairement le gain scientifique et logiciel attendu : par exemple, ameliorer la stabilite calculee, rapprocher les predictions d'un etalon, integrer un algorithme externe ou rendre la configuration plus flexible. Decrire :
- les hypotheses biologiques (temperature, ions, modifications chimiques),
- les indicateurs de qualite vises (MCC sur appariements, F1 loop detection, delta energie),
- la facon dont ces ameliorations impacteront l'interpretation fonctionnelle downstream.
Relier explicitement cet objectif au pipeline sequence -> structure -> impact afin de montrer la contribution a la feuille de route globale du projet.

## Context
Resumer les observations actuelles : sequences ou family d'ARN problematiques, limitations detectees, outils tiers utilises pour comparaison (ViennaRNA, RNAstructure, CONTRAfold, etc.). Joindre les references bibliographiques ou preprints, et indiquer :
- taille du corpus, provenances (NCBI, Rfam, bench interne),
- types de fichiers disponibles (FASTA, dot-bracket de reference, tableaux thermodynamiques),
- contraintes experimentales ou reglementaires (temperature physiologique, organisme modele).

## Proposed Change
Decrire la strategie technique :
1. Modifications prevues dans `deepstructgenomics/rna/secondary_structure.py` (energie, heuristiques, Nussinov vs autre approche).
2. Ajustements dans `deepstructgenomics/pipeline.py` pour propager options, seeds ou caches.
3. Evolution des sorties `deepstructgenomics/reporting/report_generator.py` (nouvelles colonnes, graphiques ASCII).
4. Parametres CLI supplementaires dans `scripts/run_pipeline.py` et documentation associee.
5. Eventuels ajouts dans `deepstructgenomics/config.py` pour exposer les hyperparametres.
Inclure snippets, pseudo-code et risques identifiees (complexite, memoire, precision numerique).
Preciser comment seront orchestres les benchmarks (CI, notebook partage, script shell) et quels KPI permettront d'autoriser la mise en production (variation minimale d'energie, tolerance sur le temps de calcul, etc.).

## Repro Steps / Example Inputs
Indiquer comment reproduire l'etat de l'art puis tester la proposition :
- commandes exactes (`python scripts/run_pipeline.py --sequence ...`, `--accession ...`),
- sequences ou accessions utilisees et leur longueur,
- seeds, hyperparametres, plateforme (CPU/GPU, OS, RAM),
- scripts/notebooks supplementaires si necessaire,
- methodes de comparaison (score dot-bracket, RMSD, delta energie, difference sur stability_score),
- artefacts a joindre (logs, resultats JSON, captures de rapport Markdown).
Ajouter un scenario de regression (ex : sequence connue pour echouer) et un scenario de succes attendu pour prouver la robustesse.
Si des dependances externes (ViennaRNA, outils C/C++) sont requises, preciser les commandes d'installation, les versions et les licences.

## Expected Outputs
Preciser les resultats attendus : distribution du `stability_score`, evolution du ratio de paires correctes, temps d'execution cible, niveau de confiance. Indiquer comment ces valeurs apparaitront dans les rapports JSON/Markdown et comment verifier qu'elles sont conformes (tests automatises, comparaisons externes, seuils d'acceptation). Mentionner tout changement de schema JSON (nouvelles cles, types) et si un changelog doit etre fourni aux utilisateurs. Indiquer les tickets/PR aval a reviser si la structure des sorties change (dashboards, notebooks, docs).

## Affected Modules
- `deepstructgenomics/rna/secondary_structure.py`
- `deepstructgenomics/pipeline.py`
- `deepstructgenomics/reporting/report_generator.py`
- `deepstructgenomics/config.py`
- `scripts/run_pipeline.py`
- (optionnel) dossiers `tests/`, `data/examples/`, documentation si des fichiers y sont ajoutes.

## Acceptance Criteria
- Benchmarks reproductibles fournis (script ou notebook partage dans le repo ou en annexe).
- Tests unitaires/fonctionnels crees ou ajoutes pour couvrir les nouveaux chemins.
- Execution reussie sur les sequences de demonstration existantes sans regression de temps.
- Comparaison chiffrable avec au moins un outil externe ou un dataset de reference.
- Documentation (README, docstrings) mise a jour et claire sur les nouveaux parametres.
- Rapports JSON/Markdown enrichis verifies manuellement ou automatiquement.
- Communication des resultats (issue, PR, billet interne) comprenant un resume scientifique + un resume engineering.

## Constraints / Principles
- Reproductibilite garantie (seeds fixes, parametres exposes, instructions claires).
- Pas de dependances proprietaires; privilegier bibliotheques open-source maintenables.
- Temps de calcul compatible avec un pipeline CI (indiquer les besoins si >10 min).
- Respect des choix ASCII et des conventions de logging deja en place.
- Traçabilite : decrire comment seront stockes les resultats intermediaires et comment les auditer.
- Mentionner les limites connues (pseudo noeuds, boucles multiples) et comment elles sont gerees ou ignorees.
- Specifier si des donnees sensibles ou regulatoires sont impliquees et comment elles seront protegees.
- Garder une approche modulaire afin de pouvoir remplacer le moteur de repliement sans refactor massif.

## Notes / References
Lister les DOI, liens vers preprints, billets techniques, tickets precedents. Mentionner les contributeurs/responsables a contacter, les branches existantes ou PR initiales. Ajouter, si disponible, une estimation de roadmap (jalons, dependances) et preciser ou deposer les jeux de donnees complementaires (`data/examples/`, depot externe, archive publique). Penser a inclure des captures ou extraits d'analyses qui aident a contextualiser l'issue, ainsi que les etiquettes GitHub recommandees. Si des validations externes (collaborations, comites scientifiques) sont attendues, specifier le point de contact.
Terminer par un court resume (3 lignes max) pour aider les mainteneurs a prioriser l'issue parmi les autres demandes en file d'attente.
Merci d'ajouter les labels adequats (ex: `rna`, `benchmark`, `enhancement`) lors de la soumission.
- mentionner si des jeux de donnees internes ou confidentiels existent et comment les partager sans enfreindre les licences,
- decrire les retours utilisateurs (chercheurs, medecins, bio-informaticiens) ayant motive l'issue.
