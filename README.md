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



# TODOs:
1. The database should be built from __realigned__ proteins
2. Add C2H2s!



# Owenia

```bash
bash tf_class.sh -i data/toclassify/Owefus_long.pep.fasta -o results/Owefus.db1.tsv -n 5 -d db_v1.0/
```

