for f in /oak/stanford/groups/relman/gcloud/compressed_raw_data_backup/kxue/231030-e0049-e0050-doseDependence/16S_*_R1_*.fastq.gz
do
  sample=${f##*/}
  sample=${sample%%-*}
  sample=${sample%_*}
  sample=${sample/_/-}
  sample=${sample/_/-}
  echo $f ${f/_R1_/_R2_} $sample
done > config/fastq-e0049-e0050.txt
