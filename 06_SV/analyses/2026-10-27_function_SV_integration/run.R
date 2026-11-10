#!/usr/bin/env Rscript
# Read-only integration of published-style HUMAnN pathway and SV evidence.
.libPaths(c('D:/R/library',.libPaths()))
suppressPackageStartupMessages(library(data.table))
ROOT <- 'H:/My Drive/proj_ob_AA_CRC/analyses/2026-10-27_function_SV_integration'
PROJECT <- 'H:/My Drive/proj_ob_AA_CRC'
P <- c(annotation=file.path(PROJECT,'02.processed_data/05_SV_result/SV_annotation/BMI_enriched_SV/18_SV_CDS_eggNOG_annotation_long.tsv'),
 l1=file.path(PROJECT,'analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/tables/L1_L2_candidate_master.tsv'),
 specific=file.path(PROJECT,'analyses/2026-10-13_SV_abundance_sensitivity/tables/formal_effect_specific_candidate_reporting.tsv'),
 effects=file.path(PROJECT,'analyses/2026-10-13_SV_abundance_sensitivity/tables/run_model_results_long.tsv'),
 pathways=file.path(PROJECT,'analyses/2026-06-25_high-esi-pathway-contributors/tables/high_esi_obese_crc_aligned_pathways_for_abundance_panels.tsv'),
 contributors=file.path(PROJECT,'analyses/2026-06-25_high-esi-pathway-contributors/tables/pathway_taxon_group_contribution.tsv'),
 publication=file.path(PROJECT,'analyses/2026-06-25_high-esi-pathway-contributors/tables/publication_ObAA_pathway_contributor_summary.tsv'),
 stage_means=file.path(PROJECT,'analyses/2026-09-29_BMI_enriched_SV_function_modelQC/2026-09-30_BMI_specificity_layers/data_for_figures/Heatmap_L1_dSGV_allBMI_group_means.tsv'))
stopifnot(all(file.exists(P)));hash<-tools::md5sum(P)
a<-fread(P['annotation']);l1<-fread(P['l1'])[BMI_enriched_group=='Ob']
specific<-fread(P['specific'])[BMI_enriched_group=='Ob']
effects<-fread(P['effects'])[BMI=='Ob' & Model=='abundance_adjusted']
pathways<-fread(P['pathways'])[bmi=='Ob'];pathways<-unique(pathways,by='pathway_id')
stopifnot(nrow(l1)==100L,nrow(specific)==31L,nrow(pathways)==21L)
sv<-merge(l1[,.(BMI_SV_ID,SGB,Species,SV_type,SV_feature,Original_L2_nominal=L2_nominal)],
 effects[,.(BMI_SV_ID,N,N_CTRL,N_CRC,Disease_Beta,Disease_SE,Disease_P)],by='BMI_SV_ID',all.x=TRUE)
sv[,Latest_Ob_specific:=BMI_SV_ID %in% specific$BMI_SV_ID]
# Exact SGB identity links: functional contributors are community-level pathway contributions.
c<-fread(P['contributors'])[group9=='Ob-AA']
c[,SGB:=sub('.*(SGB[0-9]+).*','t__\\1',taxon)]
c[,Pathway_ID:=sub(':.*$','',pathway_id)]
c<-c[Pathway_ID %in% pathways$pathway_id]
sv_taxa<-sv[,.(N_Ob_L1_SV=.N,N_latest_specific=sum(Latest_Ob_specific)),by=.(SGB,Species)]
linked<-merge(c,sv_taxa,by='SGB',all=FALSE)
setorder(linked,Pathway_ID,rank_abs_within_pathway)
fwrite(linked,file.path(ROOT,'tables/pathway_contributors_with_Ob_SV_evidence.tsv'),sep='\t')
# Preserve the exact overlapping CDS fraction; a partial hit is not whole-gene loss.
focus<-c('ObSV_036','ObSV_037','ObSV_038','ObSV_044','ObSV_049','ObSV_062','ObSV_130','ObSV_005')
genes<-merge(a[BMI_SV_ID %in% focus,.(BMI_SV_ID,Protein_key,CDS_ID,CDS_global_start,CDS_global_end,
 Description,Preferred_name,EC,KEGG_ko,KEGG_Pathway,PFAMs,CAZy,Max_CDS_overlap_fraction,Total_overlap_bp)],sv,by='BMI_SV_ID')
genes[,Overlap_caution:=fifelse(Max_CDS_overlap_fraction<.5,'Less than half CDS overlaps; do not infer complete gene deletion','SV overlaps CDS; annotation is predicted function, not measured activity')]
setorder(genes,BMI_SV_ID,CDS_global_start)
fwrite(genes,file.path(ROOT,'tables/focal_SV_gene_annotation_and_adjusted_effects.tsv'),sep='\t',na='NA')
# Human-audited connections: explicitly distinguish taxon overlap, function theme and exact pathway mapping.
links<-data.table(
 Topic=c('Cell surface / LPS','Species-level contributor overlap','Carbohydrate acquisition','Purine biosynthesis','Nucleotide and membrane functions','Fatty-acid biosynthesis','Histidine degradation'),
 SV_IDs=c('ObSV_037;ObSV_038','ObSV_037;ObSV_044;ObSV_049','ObSV_044;ObSV_049','ObSV_062','ObSV_036','ObSV_130;ObSV_005','ObSV_059'),
 Functional_evidence=c('Ob early shifts: PWY-7332, PWY-1269, NAGLIPASYN-PWY, PWY-8073 all increased',
 'P. vulgatus, B. stercoris and B. uniformis contribute to early-shifted pathways',
 'Bacteroides contributors recur across early-shifted pathways, including carbon-related pathways',
 'PWY-841 purine de novo biosynthesis increases in Ob AA and CRC',
 'Broad nucleotide / lipid functional remodeling; no asserted identity with PWY-841 or lipid IVA',
 'PWY66-430 myristate biosynthesis increases in Ob AA and CRC',
 'HISDEG-PWY and PWY-5030 increase in Ob AA and CRC'),
 SV_evidence=c('P. vulgatus O-antigen ligase region deletion increases; P. dorei deoxy-sugar region deletion decreases',
 'Same SGBs carry Ob SV: P. vulgatus SGB1814 specific, B. stercoris SGB1830 L1, B. uniformis SGB1836 specific',
 'B. stercoris and B. uniformis SVs overlap SusC/SusD-related proteins; deletion increases',
 'P. merdae L1 multi-segment dSV overlaps purH2/K00602 by 11.48%; deletion decreases after abundance adjustment',
 'O. splanchnicus specific vSV overlaps cls/K06131 and thyA/K00560; vSV coefficient is negative',
 'fabH/K00648-specific dSV and F. plautii L1 reductase-containing dSV both show increased deletion',
 'Histidine kinase / two-component SV annotations represent signaling, not histidine degradation'),
 Evidence_level=c('Functional-theme convergence plus some same-taxon evidence','Exact SGB overlap, not exact pathway-gene mapping',
 'Same taxa and substrate-acquisition theme, not same enzymatic pathway','Exploratory direction-compatible gene-pathway clue with weak partial overlap',
 'Specific candidate carries predicted nucleotide and phospholipid functions; direction not directly interpretable as pathway abundance',
 'Shared broad theme, but deletion direction does not support higher biosynthetic capacity','No direct SV corroboration of histidine-degradation genes'),
 Claim_boundary=c('O-antigen ligation is not lipid IVA synthesis; PWY-7332 is unclassified in contributor table; no claim of increased LPS activity',
 'Species contribution is abundance, not evidence that the SV causes its pathway signal',
 'SusC/SusD substrate cannot be inferred here; cannot infer increased carbohydrate utilization from increased deletion',
 'Not latest Ob-specific; do not call complete purH2 deletion or increased enzyme activity',
 'thyA is pyrimidine, not purine; cls is cardiolipin, not lipid IVA; vSV is not expression or activity',
 'No direct confirmation of increased myristate; host-style MetaCyc pathway names do not establish organelles in bacteria',
 'Do not equate histidine kinase / acid phosphatase with Hut-dependent histidine catabolism'))
fwrite(links,file.path(ROOT,'tables/function_SV_connection_evidence.tsv'),sep='\t')
# Check named contributor taxa missing from the Ob SV results.
named<-data.table(SGB=c('t__SGB1814','t__SGB1830','t__SGB1844','t__SGB1855'),
 Manuscript_taxon=c('Phocaeicola vulgatus','Bacteroides stercoris','Bacteroides cellulosilyticus','Bacteroides fragilis'))
named<-merge(named,sv_taxa,by='SGB',all.x=TRUE)
named[is.na(N_Ob_L1_SV),`:=`(N_Ob_L1_SV=0L,N_latest_specific=0L)]
fwrite(named,file.path(ROOT,'tables/manuscript_taxa_Ob_SV_overlap.tsv'),sep='\t')
means<-fread(P['stage_means'])[BMI_SV_ID %in% focus,.(BMI_SV_ID,Ob_CTRL,Ob_AA,Ob_CRC)]
fwrite(merge(means,sv,by='BMI_SV_ID'),file.path(ROOT,'tables/focal_dSV_stage_deletion_means.tsv'),sep='\t')
fwrite(pathways,file.path(ROOT,'tables/Ob_early_shifted_pathway_source_subset.tsv'),sep='\t')
fwrite(data.table(Input=names(P),Path=unname(P),MD5=unname(hash)),file.path(ROOT,'tables/input_hashes.tsv'),sep='\t')
stopifnot(identical(unname(hash),unname(tools::md5sum(P))))
writeLines(unname(P),file.path(ROOT,'input_files.txt'))
writeLines(c(file.path(ROOT,c('README.md','run.R','input_files.txt')),list.files(file.path(ROOT,'tables'),full.names=TRUE)),file.path(ROOT,'output_files.txt'))
print(named)
cat('Relevant PWY-841 contributors with Ob SV evidence:\n')
print(linked[Pathway_ID=='PWY-841' & SGB %in% c('t__SGB1949','t__SGB1814','t__SGB1815'),.(taxon_label,rank_abs_within_pathway,mean_pathway_fraction)])
cat('FUNCTION_SV_INTEGRATION_COMPLETE\n')
