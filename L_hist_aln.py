#!/usr/bin/env python
"""
L_hist_aln.py - library and examples to make histone sequences alignments
and to view them with annotation.

We introduce a new option here to output alignment in html.
"""
__author__="Alexey Shaytan"

import sys
sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS/')

import pandas as pd
import cPickle as pickle
from ete2 import NCBITaxa
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import json
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import uuid
from Bio.Alphabet import IUPAC

import os


from pprint import pprint


from L_aln_tools import muscle_aln, trim_aln_to_seq, trim_aln_to_seq_length
from L_shade_hist_aln import get_pdf
from hist_ss import get_hist_ss_in_aln_for_html, identify_hist_type

from L_aln2html import aln2html

os.environ['PATH']='/Users/alexeyshaytan/soft/mview-1.60.1/bin:/Users/alexeyshaytan/soft/x3dna-v2.1/bin:/Users/alexeyshaytan/soft/amber12/bin:/Users/alexeyshaytan/soft/sratoolkit/bin:/Users/alexeyshaytan/soft/bins/gromacs-4.6.3/bin:/opt/local/bin:/opt/local/sbin:/Users/alexeyshaytan/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/usr/local/ncbi/blast/bin:/usr/texbin'

# Entrez.email = "alexey.shaytan@nih.gov" 

ncbi = NCBITaxa()

def mview_plot(msa,filename):
    """ Plots an msa with mview """
    n=str(uuid.uuid4())+'.fasta'
    AlignIO.write(msa, n, "fasta")

    os.system('mview -in fasta -ruler  on -html head -coloring consensus -consensus on %s > %s'%(n,filename))
    # print os.system('echo $PATH')
    os.system("rm %s"%(n))



def annotate_hist_msa(msa,htype,variant=None):
    """Adds to the MSA lines from features.json"""

    #read json
    with open('inp_data/features.json') as ff:    
        f = json.load(ff)
    f=f[htype]
    genseq=f['General'+htype]['sequence']
    genf=f['General'+htype]['feature1']
    
    a=SummaryInfo(msa)
    cons=a.dumb_consensus(threshold=0.1, ambiguous='X')
    sr_c=SeqRecord(id='consensus',seq=cons)
    sr_genseq=SeqRecord(id='template',seq=Seq(genseq))
    auxmsa=muscle_aln([sr_c,sr_genseq])
    auxmsa.sort()

    gapped_template=str(auxmsa[1].seq)
    gapped_cons=str(auxmsa[0].seq)

    s=list()
    for c,i in zip(gapped_cons,range(len(gapped_template))):
        if(c!='-'):
            s.append(gapped_template[i])
    newgapped_template=''.join(s)
    #now we need to gap feature
    gapped_genf=list()

    k=0
    for c,i in zip(newgapped_template,range(len(newgapped_template))):
        if(c != '-'):
            gapped_genf.append(genf[i-k])
        else:
            k=k+1
            gapped_genf.append('-')
    gapped_genf=''.join(gapped_genf)


    newmsa=MultipleSeqAlignment([SeqRecord(id='gi|features|id',description=htype,seq=Seq(gapped_genf))])
    newmsa.extend(msa)
    # print newmsa
    return newmsa
    # pprint(genfeatures)



def trim_hist_aln_to_core(msa):
    """Trims hist alignment to core"""
    templ_H3 = Seq("ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVALFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA", IUPAC.protein)
    templ_H4 = Seq("SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKVFLENVIRDAVTYTEHAKRKTVTAMDVVYALKRQGRTLYGFGG", IUPAC.protein)
    templ_H2A = Seq("SGRGKQGGKTRAKAKTRSSRAGLQFPVGRVHRLLRKGNYAERVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLPKKTESSKSKSK", IUPAC.protein)
    templ_H2B = Seq("AKSAPAPKKGSKKAVTKTQKKDGKKRRKTRKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK", IUPAC.protein)

    templ_core_H3=templ_H3[43:114]
    templ_core_H4=templ_H4[23:93]
    templ_core_H2A=templ_H2A[15:119]
    templ_core_H2B=templ_H2B[33:120]

    templ={'H3':templ_core_H3,'H4':templ_core_H4,'H2A':templ_core_H2A,'H2B':templ_core_H2B}

    a=SummaryInfo(msa)
    cons=a.dumb_consensus(threshold=0.1, ambiguous='X')

    return trim_aln_to_seq_length(msa,templ[identify_hist_type(cons)])


if __name__ == '__main__':
    #Here are some examples to test this library
    #they will be also useful as examples to align and plot

    #1. Get data and prepare a subset
    #################
    hist_df=pd.read_csv('inp_data/seqs.csv') #Histone types info
    fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences
    f_hist_df=hist_df[(hist_df['hist_var']=='canonical_H2B')&(hist_df['curated']==False)]
    f_hist_df=f_hist_df.drop_duplicates(['taxid','hist_var'])[0:200]
    f_fasta_dict={key: value for (key,value) in fasta_dict.iteritems() if key in list(f_hist_df['gi'])} # get fasta dict
    # relabel with arbitrary index
    f_fasta_dict_rel={key: SeqRecord(id=str(index),seq=f_fasta_dict[key].seq) for (index,key) in enumerate(f_fasta_dict) }

    print len(f_fasta_dict)

    #2. Make MSA using my function
    #################
    # msa=muscle_aln(f_fasta_dict_rel.values()) #function takes a list of sequence records!!! #ACTIVATE FOR TEX
    msa=muscle_aln(f_fasta_dict.values()) #function takes a list of sequence records!!! #ACTIVATE FOR TEX
    AlignIO.write(msa, "int_data/msa.fasta", "fasta")
    

    #3. Get an annotated PDF of histone alignment using TEXSHADE - old way
    ##############
    #get_pdf(hist_name,align,title,shading_modes=['similar'],logo=False,hideseqs=False,splitN=20,setends=[],ruler=False):
    #The sequence names should be unique and without '|'
    if(0):
        get_pdf('H2B',msa,'H2B aln',logo=True,ruler=True)

    #4.output to html

    aln2html(msa,'int_data/h2b.html',features=get_hist_ss_in_aln_for_html(msa,'H2B',1))

    #5. TEST IT: Annotate our MSA using features.json - new experimental way
    #################
    if(0):
        annot_msa=annotate_hist_msa(msa,'H2B')
        mview_plot(annot_msa,'int_data/h2b.html')


    #6. An example of how to map alignment on one sequence
    if(0):
        templ_H2B = Seq("AKSAPAPKKGSKKAVTKTQKKDGKKRRKTRKESYAIYVYKVLKQVHPDTGISSKAMSIMNSFVNDVFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK", IUPAC.protein)
        newmsa=trim_aln_to_seq(msa,templ_H2B)
        annot_msa=annotate_hist_msa(newmsa,'H2B')
        mview_plot(annot_msa,'int_data/h2b.html')


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
    
    
            