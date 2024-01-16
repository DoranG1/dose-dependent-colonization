# Parse command-line arguments.
# -x: full path to read 1
# -y: full path to read 2
# -o: path to output directory (script creates R1 and R2 folders)
# -s: sample name (script will append -L# for each phase)
# -i: input index file (four columns: phase, read1bc, read2bc, combinedbc)
while getopts x:y:o:s:i: flag
do
	case "${flag}" in
        x) inputdemuxr1=${OPTARG};;
        y) inputdemuxr2=${OPTARG};;
        o) paramsdemuxdir=${OPTARG};;
		s) wildcardssample=${OPTARG};;
		i) inputindices=${OPTARG};;
    esac
done

#inputdemuxr1="/scratch/groups/relman/kxue/16S-pipeline/demux/16S-A1-R1.fastq.gz"
#inputdemuxr2="/scratch/groups/relman/kxue/16S-pipeline/demux/16S-A1-R2.fastq.gz"
#paramsdemuxdir="/scratch/groups/relman/kxue/16S-pipeline/demux/"
#wildcardssample="16S-A1"
#inputindices="config/indexfordemux.txt"

# Generate output directories.
mkdir -p ${paramsdemuxdir}/R1
mkdir -p ${paramsdemuxdir}/R2
# Read in the list of barcodes to parse.
while read phase r1 r2 bc
do
  # Unzip FASTQ files and extract sequencing reads with the designated barcode.
  zcat ${inputdemuxr1} | grep -A 3 _${bc} | grep -v "^--$" \
	> ${paramsdemuxdir}/R1/${wildcardssample}-L${phase}.fastq
  zcat ${inputdemuxr2} | grep -A 3 _${bc} | grep -v "^--$" \
	> ${paramsdemuxdir}/R2/${wildcardssample}-L${phase}.fastq
  # Zip output FASTQ files.
  gzip -f ${paramsdemuxdir}/R1/${wildcardssample}-L${phase}.fastq
  gzip -f ${paramsdemuxdir}/R2/${wildcardssample}-L${phase}.fastq
done < <( tail -n +2 ${inputindices} )
