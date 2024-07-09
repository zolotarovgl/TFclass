
Database v1.0

```bash
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
```

Database with realignment

__CAVE: HG realignment will take a very long time - 1000s of sequences__
Use available aligned files instead

```bash
cp ~/ant/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/data_annotation/species_list_annotation.fasta.gz .
gunzip species_list_annotation.fasta.gz
mv species_list_annotation.fasta tmp
samtools faidx tmp/species_list_annotation.fasta

# for each protein: HG.OG
cat ~/ant/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/results_annotation/gene_trees/*groups.csv  | cut -f 1-2 > tf_groups


paste <(cut -f 1 tf_groups) <(cut -f 2 tf_groups | cut -f 1 -d : | sed -E 's/\.[0-9]+$//g') <(cut -f 2 tf_groups | cut -f 1 -d :) > _tf_groups
mv _tf_groups tf_groups

N_TF=$(wc -l tf_groups  | awk '{print $1}')
N_HG=$(cut -f 2 tf_groups | sort | uniq | wc -l | awk '{print $1}')
N_OG=$(cut -f 3 tf_groups | sort | uniq | wc -l | awk '{print $1}')
echo -e "N TFs: ${N_TF}\nN HG: ${N_HG}\nN OG: ${N_OG}"

INFASTA=tmp/species_list_annotation.fasta
ALN_DIR=~/ant/gzolotarov/projects/2021_TFevol/metazoan_tf_evol_2022/030523_phylogenies/results_annotation/alignments
DBDIR=db_v2
ALN_CPU=10
mkdir -p $DBDIR
echo -e "$(date '+%Y-%m-%d %H:%M:%S')\nInput file: ${INFASTA}\nN TFs: ${N_TF}\nN HG: ${N_HG}\nN OG: ${N_OG}" > $DBDIR/info.txt

# HG level
for ID in $(cut -f 2 tf_groups  | sort | uniq); do
	echo $ID
	mkdir -p tmp/hg
	# build an hmm
	hmmbuild -n $ID tmp/hg/${ID}.hmm ${ALN_DIR}/${ID}.domains.lt.fasta > /dev/null 2>&1
done
cat tmp/hg/*.hmm > HGdb.hmm
hmmpress HGdb.hmm
mv HGdb.* $DBDIR

# OG level 
for ID in $(cut -f 3 tf_groups  | sort | uniq ); do
	echo $ID
	mkdir -p tmp/aln
	mkdir -p tmp/og
	# fetch the ids
	awk -v ID=${ID} '$3==ID{print $1}' tf_groups > tmp/aln/${ID}.ids
	N=$(wc -l tmp/aln/${ID}.ids | awk '{print $1}')
	echo "Realigning ${N} sequences"
	xargs samtools faidx ${INFASTA} < tmp/aln/${ID}.ids > tmp/aln/${ID}.fasta
	mafft --quiet --thread $ALN_CPU tmp/aln/${ID}.fasta > tmp/aln/${ID}.aln.fasta
	# build an hmm 
	hmmbuild -n $ID tmp/og/${ID}.hmm tmp/aln/${ID}.aln.fasta > /dev/null 2>&1
done
cat tmp/og/*.hmm > OGdb.hmm
#rm -rf tmp/og/
hmmpress OGdb.hmm
mv OGdb.* $DBDIR

```
