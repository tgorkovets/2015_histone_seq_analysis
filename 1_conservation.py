# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:09:07 2013

@author: alexeyshaytan

This script does conservation analysis on different subsets of histone dataset.
Outputs plots and logoPlots.

"""
from __future__ import division
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import csv
import collections
from Bio import Entrez
import cPickle as pickle
from Bio import SeqIO
  
import uuid

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess
import pylab
import pandas as pd
# import networkx as nx
from StringIO import StringIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
#from pylab import *
import itertools
import time
from multiprocessing import Pool, current_process
import cPickle as pickle

sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS/')

# from hist_ss import get_hist_ss
from L_aln_tools import muscle_aln

Entrez.email = "alexey.shaytan@nih.gov" 

tax_dict_pc={}
tax_dict_cp={}

tax_ranks_dict={}

#expand taxonomy recursive function
def get_leaves(taxids):
    """ get leaves of taxonomy given list of nodes"""
    tis=list()
    for ti in taxids:
        if ti not in tax_dict_pc:
            tis.append(ti)
        else:
            tis.extend(get_leaves(tax_dict_pc[ti]))
    return tis

def get_tree_nodes(taxids):
    """ get all nodes of taxonomy given list of nodes"""
    tis=list()
    for ti in taxids:
        if ti not in tax_dict_pc:
            tis.append(ti)
        else:
            tis.append(ti)
            tis.extend(get_tree_nodes(tax_dict_pc[ti]))
    return tis

def traverse_up(taxid):
    return [1] if taxid == 1 else [taxid]+traverse_up(tax_dict_cp[taxid])

def get_next_hub_pc(tree_pc,start,taxids):
    if start in taxids:
        return start
    if(len(tree_pc[start])==1):
        return get_next_hub_pc(tree_pc,tree_pc[start][0],taxids)
    else:
        return start

def remove_non_hubs_pc(tree_pc,start,taxids):
    d=dict()
    if(start not in tree_pc):
        return {}
    if(len(tree_pc[start])>1):
        d.update({start:tree_pc[start]})
        for i in tree_pc[start]:
            if i in tree_pc:
                d.update(remove_non_hubs_pc(tree_pc,i,taxids))
    else:
        if start in taxids:
            d.update({start:tree_pc[start]})
            # print start,' ',tree_pc[start]
            d.update(remove_non_hubs_pc(tree_pc,tree_pc[start][0],taxids))
        else:
            id=get_next_hub_pc(tree_pc,start,taxids)
            # print "slipping to ",id
            d.update({start:[id]})
            d.update(remove_non_hubs_pc(tree_pc,id,taxids))
    return d

def remove_non_hubs_cp(tree_cp_full,tree_pc,taxids):
    d=dict()
    for i in taxids:
        d.update({i:get_next_hub_cp(i,)})
    return d

def prune_tree(taxids):
    """ get clear tree for a set of taxids based on taxonomy tree"""
    #traverse up to the root - get path
    tree_pc_full=dict()
    tree_cp_full=dict()

    tree_pc=dict()
    tree_cp=dict()

    nodes=set()#get a set of all nodes that are parents of taxids plus taxids
    for ti in taxids:
        nodes=nodes.union(set(traverse_up(ti)))
    # print nodes
    #make subtree from taxonomy tree
    for i in nodes:
        if i in tax_dict_pc:
            nl=list(set(tax_dict_pc[i]).intersection(nodes))
            if len(nl):
                tree_pc_full[i]=nl
        # if in not 1:
            # tree_cp_full[i]=tax_dict_cp[i]
    # print tree_pc_full
    #remove non hub nodes, but not those in taxids
    tree_pc=remove_non_hubs_pc(tree_pc_full,1,taxids)

    #simply invert the tree
    for p in tree_pc:
        for c in tree_pc[p]:
            tree_cp.update({c:p})

    return (tree_pc, tree_cp)

def prune_subspecies(tree_pc,tree_cp,seqtaxids):

    new_taxids=list()
    for i in seqtaxids:
        if((tax_ranks_dict[i]=='subspecies') or (tax_ranks_dict[i]=='no rank')):
            if(tree_pc[tree_cp[i]][0]==i):
                new_taxids.append(i)
        else:
            new_taxids.append(i)
    return new_taxids



def main():

    #Getting data
    hist_df=pd.read_csv('inp_data/seqs.csv') #Histone types info
    fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences

    #Getting data taxonomic tree into a linked dictionary
    with open('taxdmp/nodes.dmp','r') as f:
        for line in f:
            (k1,k2,rank)=line.split('\t|\t')[0:3]
            child=int(k1)
            parent=int(k2)
            if parent==child:
                continue
            if parent in tax_dict_pc:
                tax_dict_pc[parent].append(child)
            else:
                tax_dict_pc[parent]=[child] 

            tax_dict_cp[child]=parent
            tax_ranks_dict[child]=rank

    #Getting data - taxonomics names
    df=pd.read_csv('taxdmp/names.dmp',sep='\t\|\t',converters={'type': lambda x: x[0:-2]},header=None,names=['taxid','name','name2','type'])
    df=df[df.type=="scientific name"] 
    taxid2names=pd.Series(df.name.values,index=df.taxid).to_dict() #taxid to taxname dict
    taxid2names={key:(value.split()[0][0]+'. '+' '.join(value.split()[1:2])) for (key,value) in taxid2names.iteritems()}
    

    #Here we do filtering to get a set of desired gis
    # f_hist_df=hist_df[(hist_df['curated']==True) & (hist_df['hist_type']=='H2A')]
    # f_hist_df=hist_df[(hist_df['hist_type']=='H2A')]
    f_hist_df=hist_df[(hist_df['hist_var']=='canonical_H2B')]

    #select one variant per taxid
    f_hist_df=f_hist_df.drop_duplicates(['taxid','hist_var'])


    #overlay filtering by taxids
    parent_nodes=[1] #taxids of the parent nodes we want to include.
    taxids=get_tree_nodes(parent_nodes)
    f_hist_df=f_hist_df[f_hist_df['taxid'].isin(taxids)]

    seqtaxids=list(f_hist_df['taxid'])
    print seqtaxids

    #We do not want subspecies
    #Generate a tree from available sequence taxids using original taxonomy as a guide
    #and take only one representative per species
    print "starting tree pruning"
    tree_pc,tree_cp=prune_tree(seqtaxids)
    print(tree_pc)
    print(tree_cp)
    #We need for every species only on taxid.
    print "%d sequences filtered"%len(f_hist_df)

    new_taxids=prune_subspecies(tree_pc,tree_cp,seqtaxids)
    f_hist_df=f_hist_df[f_hist_df['taxid'].isin(new_taxids)]

    print "%d sequences filtered"%len(f_hist_df)
    # exit()
    #get a list of desired fasta seqs
    f_fasta_dict={key: value for (key,value) in fasta_dict.iteritems() if key in list(f_hist_df['gi'])}
    # print f_hist_df.loc[f_hist_df.gi=='15219078','hist_type'].values[0]
    
    #Relabel sequences gi=> type and organism
    #with gi
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
    
    exit()

    #Here we construct MSA
    msa=muscle_aln(f_fasta_dict_rel.values())
    AlignIO.write(msa, "int_data/msa.fasta", "fasta")
    
    #TODO: vizualize MSA with SS annotation and conservation scores
    #TODO: map msa to reference sequence

    #We need to calculate conservation - entropy of the sequence
    #TODO: use al2co
    #TODO:
    #Here we calculate profiles


    #pylab.plot(prof)
if __name__ == '__main__':
    main()
    
            