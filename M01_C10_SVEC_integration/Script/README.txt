1.1_Transcriptomics_processing.R
Preprocessing script to prepare transcriptomic data for MOFA+ integration/SCOT+ alignment.
Input: C10SVEC_singlecells_Counts.txt
Output: processed_transcriptomics.csv, processed_transcriptomics_metadata.csv

1.2_Proteomics_processing.R
Preprocessing script to prepare proteomic data for MOFA+ integration/SCOT+ alignment.
Input: C10SVEC_singlecells_Protein_intensities.tsv, convert_uniprot.tsv
Output: processed_proteomics.csv, processed_proteomics_metadata.csv

2_SCOTplus_alignment.ipynb
SCOT+ alignment, transcriptomics projected into proteomics domain.
Input: processed_transcriptomics.csv, processed_transcriptomics_metadata.csv, processed_proteomics.csv, processed_proteomics_metadata.csv
Output: epsilon_grid_search_log.csv, alpha_comparison_summary.csv, final_foscttm_summary.csv, aligned_transcriptomics.csv, aligned_proteomics.csv, aligned_metadata.csv, sample_coupling_matrix.csv, feature_coupling_matrix.csv, Figure 1B-E, S1, S2

3_MOFA+_model_training.ipynb
MOFA+ model training script for Paired, Unpaired and SCOT+-aligned model.
Input: (Paired & Unpaired) Processed_transcriptomics.csv, Processed_proteomics.csv; (SCOT+) aligned_transcriptomics.csv, aligned_proteomics.csv
Output: (Paired) Paired_model.hdf5; (Unpaired) Unpaired_model.hdf5; (SCOT+) SCOT+-aligned_model.hdf5

4_Exploratory_downstream.R
Exploratory downstream analysis for MOFA+ integrated models.
Input: (Paired) Paired_model.hdf5; (Unpaired) Unpaired_model.hdf5; (SCOT+) SCOT+-aligned_model.hdf5
Output: Foundation figures that were used to build Figure 2-4, S3-5

5_GSEA.R
Gene set enrichment analysis for MOFA+ integrated models.
Input: (Paired) Paired_model.hdf5; (Unpaired) Unpaired_model.hdf5; (SCOT+) SCOT+-aligned_model.hdf5
Output: Foundation figures that were used to build Figure 2-4, S3-5

F2_S3_S4_S5_plotting.R, F3_plotting.R, F4_plotting.R
Plotting scripts for the corresponding figures.
Input: (Paired) Paired_model.hdf5; (Unpaired) Unpaired_model.hdf5; (SCOT+) SCOT+-aligned_model.hdf5
Output: Corresponding figures as specified by the names.

Note:
1. Skip 2_SCOTplus_alignment.ipynb for paired and unpaired workflow.
2. Figures and data are generated in the same directory as the scripts.