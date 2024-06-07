Given a protein fasta, classify the proteins into HGs and OGs. 

```bash
NCPU=30
HMMDB=db/OGdb.hmm
FASTA=data/toclassify/Octsin_long.pep.noC2H2.fasta
mkdir -p results/
OUT=results/Octsin.OG
hmmscan --cpu $NCPU -E 0.001 --tblout $OUT $HMMDB $FASTA  > /dev/null

NCPU=30
HMMDB=db/OGdb.hmm
FASTA=data/toclassify/Octsin_long.pep.noC2H2.fasta
mkdir -p results/
OUT=results/Octsin.OG
hmmscan --cpu $NCPU -E 0.001 --tblout $OUT $HMMDB $FASTA  > /dev/null
# we need to pick the best hit per query 
# awk '!x[$3]++' $OUT > $OUT.best # this one selects the best entry

```


Select the best hits to report:
```bash
Rscript select_best.r results/Octsin.HG results/Octsin.OG classification.tsv
```
# Can we add some dictionary? 
# Can we be more sure about the results? 
# We can then borrow the OG names and reference TF representatives from other species! 
# We also need to add domains! 
# We need to be more sensitive here!

# Idea 
Given a protein fasta, we need to classify the proteins into 3 levels:
1. TF Class - by PFAM: domain arrangement etc.   
2. TF Homology group - based on HMM scoring  
3. TF Orthology group - based on HMM scoring  

We can also think of an html output that will report domain this in the query and other useful information. 



