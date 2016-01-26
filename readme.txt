#Histone seq analysis project
Extract useful information form histone sequence DB.
#Data structure
inp_data - input data
L_.... - libraries with examples that perform core functions.
we might want to move some of them to MYSOFT and to CBC.
int_data - intermediate data.
taxdump - NCBI taxonomy dump.

#inp_data
seqs.csv - dump of HistoneDB 2.0 database as of October 31, 2015
with some fixes: deleted 694915110, 334187383, all NOGI
features.json - our features file from HistoneDB 2.0.

#Pipeline files
0_get_seq.py - downloads the sequences and saves them as pickle int_data/fasta_dict.p. key is the gi.


#Libraries
L_seq_subset.py - routines that subset histone sequnces, by variants, by seq features, by taxonomy.