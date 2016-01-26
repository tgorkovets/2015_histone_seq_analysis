# -*- coding: utf-8 -*-
"""
This script takes the HMM scoresa and redefines the variant, based on highest HMM.
Input: seqs_rs.cvs, scores.csv
Output: file seqs_rs_redef.csv

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
import numpy as np
  
import uuid

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess
import pylab
import pandas as pd
import networkx as nx
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
from L_aln_tools import get_prot_seqrec_by_gis

Entrez.email = "alexey.shaytan@nih.gov" 


def main():

    #1. Getting data
    df=pd.read_csv('int_data/seqs_rs.csv',dtype={'gi': 'i4'}) #Histone types info
    scores=pd.read_csv('inp_data/scores.csv',dtype={'gi': 'i4'})
    # fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences
    new_var_list=list()
    gis=list(df['gi'])
    n=len(gis)
    i=0
    for gi in gis:
        i+=1
        if(i%100==0):
            print i*100/n
        s=scores[(scores['gi']==int(gi))]
        if(len(s)>0):
            s2=s.sort(['score'],ascending=False)
            # print s2
            new_var=s2['hist_var'].values[0]
            # print new_var
        else:
            new_var=df[(df['gi']==gi)]['hist_var'].values[0]
        # print new_var
        new_var_list.append(new_var)
    df['hist_var']=new_var_list

    df.to_csv('int_data/seqs_rs_redef.csv',index=False)

if __name__ == '__main__':
    main()
    
            