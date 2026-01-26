# Validation des recuperations NCBI et enrichissement metadata

## Perimetre
- [ ] Robustesse des appels `efetch` / `esummary` dans `deepstructgenomics/data_sources/ncbi_client.py`
- [ ] Propagation des metadonnees dans `deepstructgenomics/pipeline.py`
- [ ] Formatage des rapports (`deepstructgenomics/reporting/report_generator.py`)
- [ ] Experience CLI (`scripts/run_pipeline.py`) et messages d'erreur
- [ ] Gestion du cache / dossier `data/cache/`
Choisissez les cases pertinentes et detaillez si des endpoints specifiques (nuccore, gene, snp) sont vises.
Indiquez egalement si l'issue concerne des executions planifiees (batch) ou des usages ponctuels afin de calibrer la priorisation.

## Objective
Expliciter l'objectif principal : fiabiliser les telechargements, enrichir les metadonnees exposees aux modules aval, gerer les quotas ou fournir des validations automatiques. Decrire les raisons scientifiques (annotations correctes pour l'analyse structurale, provenance tracee, conformite clinique) et les objectifs techniques (temps de reponse, tolerance aux erreurs, reessais).
Relier l'objectif aux besoins des autres modules (RNA folding, variants, reporting) afin de souligner l'impact fonctionnel.
Indiquer les livrables attendus (PR, notebook, dashboard) pour cadrer la definition du done.

## Context
Lister les accessions et bases reclamees (nuccore, nucleotide, gene). Preciser les problemes rencontrees (timeouts, HTTP 429, champs manquants, differences selon l'UID). Indiquer :
- environnement d'execution (proxy, VPN, conteneur, cluster),
- presence ou non d'une cle API et valeur de `tool`/`email`,
- historique d'incidents ou issues existantes,
- besoins specifiques en metadonnees (taxid, organism, keywords, bioproject, references).
Joindre si possible des extraits anonymises de reponses conseillées et de reponses problematiques.

## Proposed Change
Decrire les modifications envisagees :
1. Ajout de validations schema / parsing resilient dans `deepstructgenomics/data_sources/ncbi_client.py` (retry, backoff, log structure).
2. Mise a jour de `deepstructgenomics/pipeline.py` pour stocker les nouveaux champs et les rendre disponibles aux modules RNA/variants.
3. Ajustements de `deepstructgenomics/reporting/report_generator.py` afin d'inclure les champs dans le JSON et la section "Metadonnees" du Markdown.
4. Ameliorations CLI `scripts/run_pipeline.py` (options `--ncbi-email`, `--ncbi-api-key`, messages utilisateur).
5. Documentation du dossier `data/cache/` (strategie purge, TTL, taille max).
Inclure pseudo-code pour le retry/backoff, tableaux des champs cibles et risques courses (quota, latence, compatibilite offline).
Mentionner comment seront ajoutes les tests automatiques (mocks HTTP, enregistrements VCR, caches pre-remplis) et comment seront suivis les indicateurs de qualite.
Si une evolution des dependances (requests, retrypolicies) est necessaire, lister clairement les impacts packaging (requirements, pyproject, CI).

## Repro Steps / Example Inputs
Fournir un protocole reproductible :
- commandes CLI exactes (`python scripts/run_pipeline.py --accession NR_000027.1 --ncbi-email ...`),
- parametres annexes (API key, proxies, headers),
- sequences/accessions problematiques avec copies des reponses NCBI anonymisees,
- scripts de test (unitaires ou smoke) qui valident les nouveaux comportements,
- methodes pour simuler des erreurs (throttling, reponse vide, JSON invalide),
- methodes de verification (diff JSON, inspection Markdown, logs).
Ajouter un mode d'emploi pour reprendre les tests en cas de limites de quotas (delais d'attente, parametre de throttle) et un scenario offline de secours.
Preciser les commandes ou scripts pour nettoyer `data/cache/` avant et apres la reproduction.

## Expected Outputs
Decrire les nouveaux champs attendus dans les rapports (ex: `taxid`, `subtype`, `bioproject`, `keywords`, `gb_source`). Expliquer les formats (string, int, liste) et comment ils se retrouveront dans les sections JSON et Markdown. Préciser les seuils de succes : par exemple, 100 % des accessions testees contiennent un champ `organism`, un retry automatique resout >90 % des erreurs 429, le cache reduit les appels de 40 %, etc.
Indiquer si un changelog ou une migration sera requis pour les utilisateurs downstream (dashboards, outils internes) et comment la compatibilite ascendante sera preservee.
Decrire la facon dont les nouveaux champs seront valides (schéma JSON, tests pytest, diff Markdown signe).

## Affected Modules
- `deepstructgenomics/data_sources/ncbi_client.py`
- `deepstructgenomics/pipeline.py`
- `deepstructgenomics/reporting/report_generator.py`
- `scripts/run_pipeline.py`
- (optionnel) scripts/tests relies au cache ou a l'import de donnees.

## Acceptance Criteria
- Suite de tests couvrant les cas reussis, les erreurs HTTP, les reponses partielles et les contenus non ASCII.
- Documentation mise a jour (README, docstrings, gabarits CLI) pour expliquer email/API key, limites de quotas et comportement du cache.
- Logs clairs lorsqu'une requete echoue (code, retries, suggestion pour l'utilisateur).
- Demonstation de plusieurs accessions representatives (ARN ribosomique, mRNA, lncRNA) avec sortie JSON/Markdown verifiee.
- Evidence que `deepstructgenomics/pipeline.py` et `deepstructgenomics/reporting/report_generator.py` consomment bien les nouveaux champs.
- Rapport de validation resument les tests executes, leurs resultats et les eventuels suivis.

## Constraints / Principles
- Respect strict des Conditions d'utilisation NCBI : throttle, identification du tool, pas de scraping massif.
- Eviter de stocker des donnees sensibles dans `data/cache/`; decrire la strategie de purge.
- Garantir un mode degrade (message clair) lorsque le reseau ou l'API ne sont pas accessibles.
- Rester compatible avec l'execution offline si la sequence est fournie localement (ne pas bloquer le pipeline).
- Utiliser uniquement des dependances open-source et minimales pour ne pas alourdir l'installation.
- Documenter toute decision de caching ou de stockage (taille max, encryption, rotation).
- S'assurer que les messages d'erreur restent en francais/anglais clair et incluent les identifiants utiles (accession, code HTTP).

## Notes / References
Ajouter les liens vers la documentation E-utilities (Esummary/Efetch/Einfo), tickets precedents et exemples de reponses JSON. Mentionner les personnes de reference (contact NCBI, DS/IT interne). Indiquer ou stocker les traces (logs, captures) et, si pertinent, un calendrier de mise en production. Terminer par un resume synthetique pour contextualiser la priorite de l'issue. Preciser aussi les etiquettes GitHub recommandees (`ncbi`, `data-source`, `bug`, `enhancement`) et les dependances vers d'autres issues si elles existent.
Si une review legale/compliance est necessaire (gestion emails, quotas), mentionnez la personne a solliciter et l'estimation de delai.
Pensez a ajouter les labels adequats et a rappeler si l'issue est bloquante pour une autre iteration du pipeline.
Merci de joindre un court resume final pour aider a prioriser la demande.
