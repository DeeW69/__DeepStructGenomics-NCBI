# Ameliorations de l'analyse d'impact des variants

## Perimetre
- [ ] Nouvelles heuristiques ou scores dans `deepstructgenomics/variants/variant_analysis.py`
- [ ] Interaction avec la structure (`deepstructgenomics/rna/secondary_structure.py`)
- [ ] Mise a jour du pipeline (`deepstructgenomics/pipeline.py`)
- [ ] Sorties JSON/Markdown (`deepstructgenomics/reporting/report_generator.py`)
- [ ] Options CLI / entrées fichiers (`scripts/run_pipeline.py`)
Indiquez si vous traitez des SNP isoles, des indels, des variants multiples ou des fichiers externes (VCF, JSON).

## Objective
Preciser la question scientifique : detection de variants destabilisant la structure ARN, integration d'informations ClinVar, scoring quantitatif, priorisation pour analyses fonctionnelles. Mentionner :
- les types de variants (substitutions simples, indels, variations multiples),
- les hypotheses biologiques (isoformes, contexte cellulaire, temperature),
- les indicateurs de succes (delta energie, classification pathogene vs benign, gain de precision).
Relier la demande au pipeline sequence -> structure -> impact pour montrer la valeur medicinale/biomedicale.
Indiquer clairement quels acteurs (chercheurs, cliniciens, data scientists) utiliseront le nouveau score et sous quelle forme.

## Context
Donner les jeux de donnees et references :
- sources (ClinVar, dbSNP, VCF internes),
- sequences ou accessions NCBI cibles,
- outils tiers utilises pour comparaison (RNAsnp, RNAmute, autres frameworks),
- limitations actuelles observees dans `variant_analysis` (score trop simple, manque de support pour indels, etc.).
Indiquer les contraintes experimentales (longueur max, type d'organisme, besoin d'anonymisation des variants).
Si des donnees proprietaires sont impliquees, decrire comment partager des exemples synthetiques.

## Proposed Change
Decrire le plan technique :
1. Evolutions proposes dans `deepstructgenomics/variants/variant_analysis.py` (nouvelles structures de donnees, matrices energie, per-variant metadata).
2. Interaction avec `deepstructgenomics/rna/secondary_structure.py` pour recalcul partiel ou alignement structurel.
3. Adaptations de `deepstructgenomics/pipeline.py` pour accepter des fichiers variants, orchestrer plusieurs comparaisons, gerer les logs.
4. Mise a jour de `deepstructgenomics/reporting/report_generator.py` pour afficher les nouveaux champs (scores, commentaires, tableaux).
5. Ajout de parametres dans `scripts/run_pipeline.py` (chemin VCF, format JSON, options de filtrage).
Inclure pseudo-code, diagrammes ou exemples de structures, ainsi que les risques identifies.
Preciser la strategie de tests et l'impact sur les dependances (nouveaux packages, scripts).

## Repro Steps / Example Inputs
Expliquer comment reproduire l'etat actuel et valider la proposition :
- sequences reference/mutant (CLI `--sequence` ou `--accession` + `--mutant-sequence`),
- exemples de variants (HGVS, VCF, JSON) et formats d'entree,
- commandes CLI et scripts utilises,
- seeds ou parametres specifiques,
- methodes de comparaison (alignement, evaluation sur dot-bracket, verif ClinVar),
- jeux de tests existants a reutiliser.
Ajouter un scenario ou l'algorithme actuel echoue et un scenario cible brute pour quantifier le gain.
Si des pipelines externes (VCF parsing, annotation ClinVar) sont necessaires, detaillez leur version et leur licence.

## Expected Outputs
Decrire les nouveaux champs attends dans `variant_result` (scores multiples, classification, listes de paires gagnees/perdues, commentaires). Expliquer comment ces donnees apparaitront :
- dans le JSON (schema, types, exemple de payload),
- dans le Markdown (nouveaux tableaux, sections additionnelles).
Definir les seuils de succes (precision minimale, recall, temps d'execution) et comment ils seront mesures (tests automatises, evaluation manuelle).
Mentionner les besoins de migration / compatibilite (par ex. nouveaux champs optionnels vs obligation de reprocesser d'anciens rapports).

## Affected Modules
- `deepstructgenomics/variants/variant_analysis.py`
- `deepstructgenomics/rna/secondary_structure.py`
- `deepstructgenomics/pipeline.py`
- `deepstructgenomics/reporting/report_generator.py`
- `scripts/run_pipeline.py`
- (optionnel) nouveaux fichiers de configuration, jeux de donnees, tests.

## Acceptance Criteria
- Cas de test couvrant au moins un variant benign, un variant deletere et un variant incertain.
- Documentation et docstrings mises a jour (README, commentaires).
- Rapport CLI de demonstration prouvant que les nouveaux champs apparaissent correctement.
- Comparaison chiffrable avec un outil ou benchmark externe.
- Maintien de la compatibilite retro (anciens arguments/outputs fonctionnent toujours).
- Instructions claires pour importer des jeux de variants plus larges (VCF, JSON).
- Communication finale (issue, PR) contenant un resume scientifique + technique.

## Constraints / Principles
- Transparence des scores : decrire formules, hypotheses thermodynamiques ou statistiques.
- Utilisation d'outils open-source, pas de dependances proprietaires.
- Temps de calcul raisonnable (indiquer si un traitement batch ou GPU est necessaire).
- Reproductibilite (seeds, parametres exposes, instructions CLI).
- Traitement responsable des donnees patients : anonymisation, stockage temporaire, pas de fuite en logs.
- Modularite pour permettre l'ajout futur d'autres modeles (ML, deep learning) sans refactor massif.

## Notes / References
Lister les publications, DOI, billets techniques, issues liees. Mentionner les contributeurs cle et les reviewers cibles (experts variants, experts structure). Indiquer ou stocker les donnees d'exemple (dossier `data/examples/`, lien externe). Ajouter, si disponible, une roadmap (jalons, dependances) et les etiquettes GitHub souhaitees (`variants`, `impact-analysis`, `enhancement`). Terminer par un court resume pour aider a prioriser l'issue.
Ajouter les labels appropries et preciser si l'issue est bloquante pour une version ou une etape clinique.
