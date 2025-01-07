nextflow run main.nf --matrix dataset/sub_pbmc_matrice.csv \
                     --transcriptome dataset/Homo_sapiens.GRCh38.cdna.all.fa \
                     --features gene_name \
                     --sim_celltypes true \
                     --cell_types_annotation dataset/sub_pbmc_cell_type.csv \
                     --gtf dataset/GRCh38-2020-A-genes.gtf \
                     --fastq_model dataset/SC3pv3_GEX_Human_PBMC_ONT_1M.fastq \
                     --full_length false \
                     --pcr_cycles 2 \
                     --pcr_total_reads 20000000 \
                     --threads 30 \
                     --build_model true \
                     --projetctName "paper_PBMC_5k_30M"
         

