rawdir="/oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/231030-e0049-e0050-doseDependence/"
outdir="processStrainIsolates"

run_midas.py species ${outdir}/SKX01232 \
   -1 ${rawdir}/SKX01232_TAGAATTGGA-GTTAATCTGA_L001_R1_001.fastq.gz -2 ${rawdir}/SKX01232_TAGAATTGGA-GTTAATCTGA_L001_R2_001.fastq.gz -t 16
   
run_midas.py species ${outdir}/SKX01292 \
   -1 ${rawdir}/SKX01292_GATATTGTGT-ACCACACGGT_L001_R1_001.fastq.gz -2 ${rawdir}/SKX01292_GATATTGTGT-ACCACACGGT_L001_R2_001.fastq.gz -t 16
   
run_midas.py species ${outdir}/SKX01304 \
   -1 ${rawdir}/SKX01304_CACTCAATTC-TTAAGTTGTG_L001_R1_001.fastq.gz -2 ${rawdir}/SKX01304_CACTCAATTC-TTAAGTTGTG_L001_R2_001.fastq.gz -t 16
   
run_midas.py species ${outdir}/SKX01305 \
   -1 ${rawdir}/SKX01305_TTCTGTAGAA-GAGCGCAATA_L001_R1_001.fastq.gz -2 ${rawdir}/SKX01305_TTCTGTAGAA-GAGCGCAATA_L001_R2_001.fastq.gz -t 16
