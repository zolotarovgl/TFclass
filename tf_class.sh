
#! /usr/bin/bash
set -euo pipefail
FASTA=$1
OUTFILE=$2


NCPU=30
mkdir -p tmp/

echo "Running classifications using ${NCPU} threads..."

HMMDB=db/HGdb.hmm
HG_OUT=tmp/pred.HG
hmmscan --cpu $NCPU -E 0.001 --tblout $HG_OUT $HMMDB $FASTA  > /dev/null
echo "Done HG-level classifications"
HMMDB=db/OGdb.hmm
OG_OUT=tmp/pred.OG
hmmscan --cpu $NCPU -E 0.001 --tblout $OG_OUT $HMMDB $FASTA  > /dev/null
echo "Done OG-level classifications"

# Output file
Rscript select_best.r ${HG_OUT} ${OG_OUT} ${OUTFILE}

