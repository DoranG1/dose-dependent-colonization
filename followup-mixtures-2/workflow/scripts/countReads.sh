rm -f ../out/trimmed/e0049-e0050/summary.txt
for f in ../out/trimmed/e0049-e0050/R1/*.fastq.gz
do
  numReads=$( zcat $f | grep ^@ | wc -l )
  echo $f $numReads
done > ../out/trimmed/e0049-e0050/summary.txt
