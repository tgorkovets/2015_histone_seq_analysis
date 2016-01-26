# -*- coding: utf-8 -*-
"""
This script analyzes the downloaded sequences, determines if a sequence is refseq and adds this info to the table.
New file seqs_rs.csv

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
    df=pd.read_csv('inp_data/seqs.csv') #Histone types info
    fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences
    # exit()
    rs_list=list()
    p_list=list()
    for i,row in df.iterrows():
        gi=row['gi']
        id=fasta_dict[str(gi)].id
        desc=fasta_dict[str(gi)].description

        if(re.search('(XP_)',id)):
            rs=1
        else:
            if(re.search('(NP_)',id)):
                rs=2
            else:
                rs=0
        if(re.search('partial',desc)):
            p=True
        else:
            p=False
        rs_list.append(rs)
        p_list.append(p)
    df['RefSeq']=rs_list 
    df['partial']=p_list

    df.to_csv('int_data/seqs_rs.csv',index=False)

if __name__ == '__main__':
    main()
    
            