#!/usr/bin/env bash
set -euo pipefail

# Usage function
usage() {
    echo "Usage: $0 -i <input_fasta> -o <output_file> -n <num_cpu> -d <db_dir>"
    exit 1
}

# Parse command-line arguments
while getopts ":i:o:n:d:" opt; do
  case ${opt} in
    i )
      FASTA=$OPTARG
      ;;
    o )
      OUTFILE=$OPTARG
      ;;
    n )
      NCPU=$OPTARG
      ;;
    d )
      DB_DIR=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      usage
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      usage
      ;;
  esac
done
shift $((OPTIND -1))

# Check if mandatory arguments are set
if [ -z "${FASTA:-}" ] || [ -z "${OUTFILE:-}" ] || [ -z "${NCPU:-}" ] || [ -z "${DB_DIR:-}" ]; then
    echo "Error: Input fasta file, output file, number of CPUs, and database directory are mandatory"
    usage
fi

# Create temporary directory
mkdir -p tmp/

echo "Running classifications using ${NCPU} threads..."

# Run HG-level classifications
HGDB=${DB_DIR}/HGdb.hmm
HG_OUT=tmp/pred.HG

N=$(grep -c '>' $FASTA)

echo "hmmscan: ${N} proteins ..."

hmmscan --cpu $NCPU -E 0.001 --tblout $HG_OUT $HGDB $FASTA > /dev/null
echo "Done HG-level classifications"

# Run OG-level classifications
OGDB=${DB_DIR}/OGdb.hmm
OG_OUT=tmp/pred.OG
echo "hmmscan: ${N} proteins ..."
hmmscan --cpu $NCPU -E 0.001 --tblout $OG_OUT $OGDB $FASTA > /dev/null
echo "Done OG-level classifications"

# Output file
echo "Rscript select_best.r ${HG_OUT} ${OG_OUT} ${OUTFILE}"
Rscript select_best.r ${HG_OUT} ${OG_OUT} ${OUTFILE}

