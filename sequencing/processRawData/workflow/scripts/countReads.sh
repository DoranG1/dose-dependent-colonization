rm -f ../out/trimmed/e0043/summary.txt
for f in ../out/trimmed/e0043/R1/*.fastq.gz
do
  numReads=$( zcat $f | grep ^@ | wc -l )
  echo $f $numReads
done > ../out/trimmed/e0043/summary.txt

rm -f ../out/trimmed/e0063/summary.txt
for f in ../out/trimmed/e0063/R1/*.fastq.gz
do
  numReads=$( zcat $f | grep ^@ | wc -l )
  echo $f $numReads
done > ../out/trimmed/e0063/summary.txt

rm -f ../out/trimmed/e0065/summary.txt
for f in ../out/trimmed/e0065/R1/*.fastq.gz
do
  numReads=$( zcat $f | grep ^@ | wc -l )
  echo $f $numReads
done > ../out/trimmed/e0065/summary.txt

rm -f ../out/trimmed/blank/summary.txt
for f in ../out/trimmed/blank/R1/*.fastq.gz
do
  numReads=$( zcat $f | grep ^@ | wc -l )
  echo $f $numReads
done > ../out/trimmed/blank/summary.txt
