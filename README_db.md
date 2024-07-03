# HG level 
# build a HG database 
cp ~/ant/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/results_annotation/alignments/tfs.*lt.fasta data/alignments/
mkdir -p db
for FILE in $(ls data/alignments/tfs.*lt.fasta); do
echo $FILE
NAME=$(basename $FILE  | cut -f 2-3 -d .)
echo $NAME
hmmbuild -n $NAME db/${NAME}.hmm $FILE
done
cat db/*hmm  > HGdb.hmm
rm db/*hmm
hmmpress HGdb.hmm
mv HGdb.* db

# OG database
# Need a classification of available prots into OGs. 
cat ~/ant/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/results_annotation/gene_trees/tfs.*groups.csv  | cut -f 1-2 | cut -f 1 -d :  | sed 's/tfs.//g' > data/classifications/tf2og

mkdir -p tmp/OG_alignments
#OG=Homeodomains.HG9.11
for OG in $(cat data/classifications/tf2og  | cut -f 2  | sort | uniq ); do
HG=$(echo $OG | cut -f 1-2 -d .)
echo $OG
echo $HG
samtools faidx data/alignments/tfs.${HG}.domains.lt.fasta
xargs samtools faidx data/alignments/tfs.${HG}.domains.lt.fasta < <(grep -w $OG data/classifications/tf2og  | cut -f 1) > tmp/OG_alignments/${OG}.fasta
hmmbuild -n $OG tmp/OG_alignments/${OG}.hmm tmp/OG_alignments/${OG}.fasta
done
cat tmp/OG_alignments/*hmm > tmp/OGdb.hmm
hmmpress tmp/OGdb.hmm
mv tmp/OGdb.hmm* db/
#rm -rf tmp/OG_alignments/

