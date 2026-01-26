# Gabarits d'issues DeepStructGenomics-NCBI

Ce dossier regroupe les modeles Markdown utilises pour cadrer les contributions autour du pipeline sequence -> structure -> impact. Ils renvoient directement aux modules clefs (`deepstructgenomics/rna/secondary_structure.py`, `deepstructgenomics/data_sources/ncbi_client.py`, `deepstructgenomics/variants/variant_analysis.py`, `deepstructgenomics/pipeline.py`, `deepstructgenomics/reporting/report_generator.py`, `scripts/run_pipeline.py`).

## Choisir le bon modele

- `rna_folding_improvements.md` : ameliorations algorithmiques, benchmarks ou nouvelles metriques de repliement ARN.
- `ncbi_fetch_validation.md` : validation des appels NCBI, enrichissement metadata, resiliences reseau/cache.
- `variant_impact_improvements.md` : evolution des scores reliant mutation -> structure -> impact fonctionnel.

## Informations minimales a fournir

1. Objectif scientifique (hypotheses, valeur attendue) et livrables techniques.
2. Contexte / references (publications, jeux de donnees, incidents rencontres).
3. Description precise des changements proposes (modules, scripts, configurations).
4. Donnees de reproduction (commandes, sequences, variants, parametres).
5. Critere d'acceptation et contraintes (reproductibilite, performance, licences, open-source).

## Exemples de titres d'issue

- Repliement : "Benchmark Nussinov vs outil externe sur NR_000027"
- NCBI : "Audit metadonnees Entrez pour NR_003286"
- Variant : "Score delta energie pour SNP rsXXXX dans `variant_analysis`"

Selectionnez le modele adequat lors de la creation d'une issue GitHub, remplissez chaque section de maniere concise mais exploitable et joignez les ressources pertinentes (liens, jeux de tests, contacts). Cela facilite la revue scientifique comme logicielle et aligne les contributions sur la feuille de route.
