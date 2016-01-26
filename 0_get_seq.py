# -*- coding: utf-8 -*-
"""
This script just download all the sequences from our seqs.cvs file - the histone database.

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

    hist_df=pd.read_csv('inp_data/seqs.csv')
    gis=list(hist_df['gi'])
    gis=[i for i in gis if re.match('\d+',i)]
    fasta_dict=get_prot_seqrec_by_gis(gis)
    pickle.dump( fasta_dict, open("int_data/fasta_dict.p", "wb" ) )
        #print manual
        # manual['GI']=manual['GI'].astype(str)
        # human_hist_df=manual[manual['Taxa']=='human']
        # print human_hist

    #pylab.plot(prof)
if __name__ == '__main__':
    main()
    
            