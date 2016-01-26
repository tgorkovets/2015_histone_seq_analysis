#!/usr/bin/env python
"""
L_seq_subset.py - library and examples to subset histone sequences list.

Various functions to do subsetting by taxonomy, taxonomy subsampling,
filtering by sequence quality,
as well as examples on how to filter by variant.
"""
__author__="Alexey Shaytan"

import sys
import numpy as np
import pandas as pd
import cPickle as pickle
from ete2 import NCBITaxa
from pprint import pprint
import os.path
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS/')
# Entrez.email = "alexey.shaytan@nih.gov" 
from hist_ss import get_core_lendiff
rank_dict={'superkingdom':0,'kingdom':1,'phylum':2,'class':3,'superorder':4,'order':5,'suborder':6,'infraorder':7,'parvorder':8,'superfamily':9,'family':10,'subfamily':11,'genus':12,'subgenus':13,'species':14,'subspecies':15}

ncbi = NCBITaxa()

def check_hist_length(seq,hist_type,hist_var=None,dev_percent=10):
    """
    This simple check compares the length of sequence provided
    to a range of curated sequences in histone DB +- dev_percent %.
    """
    if(os.path.isfile('int_data/cur_length.csv')):
        cur_df=pd.read_csv('int_data/cur_length.csv')
    else:
        hist_df=pd.read_csv('inp_data/seqs.csv') #Histone types info
        fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences
        #construct df with length
        cur_df=hist_df[(hist_df['curated']==True)]
        cur_df['length']=cur_df['gi'].map(lambda x: len(fasta_dict[str(x)].seq))
        # print cur_df.groupby(['hist_type','hist_var']).agg([np.max,np.min])
        cur_df.to_csv('int_data/cur_length.csv')

    length=len(seq)
    if(hist_var):
        f_df=cur_df[(cur_df['hist_var']==hist_var)]
    else:
        f_df=cur_df[(cur_df['hist_type']==hist_type)]
    min_l=min(f_df['length'])
    max_l=max(f_df['length'])
    if (length>min_l*(1.0-dev_percent/100.)) and (length<max_l*(1.0+dev_percent/100.)):
        return True
    else:
        return False


def check_hist_core_length(seq,hist_type,dev_percent=5):
    """
    This simple check compares the length of core part of histone +- dev_percent %.
    """
    #Let's define 1kx5 sequences
    templ_H3 = Seq("ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVALFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA", IUPAC.protein)
    templ_H4 = Seq("SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", IUPAC.protein)
    templ_H2A = Seq("SGRGKQGGKTRAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPKKTESSKSKSK", IUPAC.protein)
    templ_H2B = Seq("AKSAPAPKKGSKKAVTKTQKKDGKKRRKTRKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK", IUPAC.protein)
    temp_seq={'H2A':templ_H2A,'H2B':templ_H2B,'H3':templ_H3,'H4':templ_H4}
    ratio=get_core_lendiff(seq,temp_seq[hist_type],hist_type,debug=0)
    # print ratio,(1.0-dev_percent/100.)
    if (ratio>(1.0-dev_percent/100.)) and (ratio<(1.0+dev_percent/100.)):
        return True
    else:
        return False
 

# def subsample_taxids(taxids,rank='species'):
#     """
#     For a given set of taxids leaves only one representative per selected rank
#     Eg. for a set of subspecies - leave only species.
#     """

#     tree = ncbi.get_topology(taxids,intermediate_nodes=True)

#     # print tree.get_ascii(attributes=["sci_name", "rank","taxid"])
#     #We have now a phylogenetic tree with all annotations for our taxids.
#     subsampled_taxids=set()
#     #We are iterating through the taxids and for every we are determining only one representative from this group.
#     #These representatives will be the same for the taxids in one group - and hence subsampling will happen.
#     for t in taxids:
#         # print "Doing for",t
#         #From a taxid we need to go up the tree till we reach the desired rank.
#         node=tree.search_nodes(name=str(t))[0]

#         # print rank_dict.get(node.rank,100),'vs',rank_dict.get(rank)

#         while rank_dict.get(node.rank,100)>rank_dict.get(rank):
#             node=node.up

#         # print "Uppernode"
#         # print node.name

#         #And now we go down taking the first child taxid.
#         while str(node.name) not in map(str,taxids):
#             node=node.children[0]
#         subsampled_taxids.add(int(node.name))
#     return list(subsampled_taxids)



if __name__ == '__main__':
    #Here are some examples to test this library
    #they will be also useful as examples to filter histone dataframe

    #1. Getting data
    #################
    hist_df=pd.read_csv('inp_data/seqs.csv') #Histone types info
    fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences

    #2. Filter dataframe by histone variant
    #################
    f_hist_df=hist_df[(hist_df['hist_var']=='canonical_H2B')]

    #3. Select one variant per taxid
    #################
    f_hist_df=f_hist_df.drop_duplicates(['taxid','hist_var'])
    
    #4. Filter by list of taxonomy clades   
    ################
    parent_nodes=[9443] #131567 - cellular organisms
    taxids=list()
    for i in parent_nodes:
        taxids.extend(ncbi.get_descendant_taxa(i))

    f_hist_df=f_hist_df[f_hist_df['taxid'].isin(taxids)]

    #5. Take one representative per species or specific rank.
    ################
    #Common ranks: superorder-order-suborder-infraorder-parvorder-superfamily-family-subfamily-genus-species-subspecies
    seqtaxids=list(f_hist_df['taxid']) #old list
    new_seqtaxids=subsample_taxids(seqtaxids,rank='family') #new subsampled list
    f_hist_df=f_hist_df[f_hist_df['taxid'].isin(new_seqtaxids)] #remake the dataframe
    
    #---------------
    #Output tree before subsampline
    tree = ncbi.get_topology(seqtaxids,intermediate_nodes=False)
    print tree.get_ascii(attributes=["sci_name", "rank","taxid"])
    #Output after subsampling
    tree = ncbi.get_topology(new_seqtaxids,intermediate_nodes=False)
    print tree.get_ascii(attributes=["sci_name", "rank","taxid"])
    #Output new and old sets
    print seqtaxids
    print new_seqtaxids

    #6. Filter sequences by their quality: align to known histones and remove those that have strange indels.
    ################
    #### it turned out many sequences are annotated by pipeline incorrectly, incorrect start and stop codons
    #### often contains two histones from cluster that are accidentally merged
    #### one way is to restrict sequences by length.
    #### the other is to restrict by core length
    
    newgis=list()
    for i,row in f_hist_df.iterrows():
        gi=row['gi']
        seq=fasta_dict[str(gi)].seq
        hist_type=row['hist_type']
        hist_var=row['hist_var']
        print check_hist_length(seq,hist_type,hist_var,10)
        print check_hist_core_length(seq,hist_type,10)

        # if(check_hist_length(seq,hist_type,hist_var,10)&check_hist_core_length(seq,hist_type,hist_var,10)):
            # newgis.append(gi)

    exit()


    #############################################################################
    #############################################################################
    #############################################################################
    #############################################################################

    #6. Tricks with data frames and dictionaries - NOT tested
    ########################

    #Relabel sequences gi=> type and organism
    f_fasta_dict_rel={key: SeqRecord(id=key, description=f_hist_df.loc[f_hist_df.gi==key,'hist_var'].values[0]+' '+taxid2names[f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]],seq=value.seq) for (key,value) in f_fasta_dict.iteritems() }
    #with arbitrary index
    # f_fasta_dict_rel={key: SeqRecord(id=str(index), description=f_hist_df.loc[f_hist_df.gi==key,'hist_var'].values[0]+' '+taxid2names[f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]],seq=f_fasta_dict[key].seq) for (index,key) in enumerate(f_fasta_dict) }
    keys=str()

    #output taxids  
    for (key,value) in f_fasta_dict.iteritems():
        keys=keys+str(f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0])+','
    print keys

    #output patternmatch H2A.Z
    # for (key,value) in f_fasta_dict.iteritems():
        # if(re.search('R[VI][GSA][ASG]K[SA][AGS]',str(value.seq))):
            # print "%s,%s"%(str(f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]),'#00ff00')
        # else:
            # if(re.search('R[VI][GSA][ASG]G[SA]P',str(value.seq))):
                # print "%s,%s"%(str(f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]),'#0000ff')
            # else:
                # print "%s,%s"%(str(f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]),'#ff0000')
    
    #output patternmatch H2B
    for (key,value) in f_fasta_dict.iteritems():
        if(re.search('[^K]$',str(value.seq))):
            print "%s,%s"%(str(f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]),'#00ff00')
        else:
            print "%s,%s"%(str(f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]),'#ff0000')
    
    #Here we construct MSA
    msa=muscle_aln(f_fasta_dict_rel.values())
    AlignIO.write(msa, "int_data/msa.fasta", "fasta")
    
    
            