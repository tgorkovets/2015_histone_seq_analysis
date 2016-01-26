# -*- coding: utf-8 -*-
"""
This script interactively explores H2A alignments and
conservation using tools and examples developed in libraries.
Use it as a starting example, to further modify it.
"""
import sys

sys.path.append('/Volumes/MDBD/Dropbox/work/MYSOFT/ALIGNMENT_TOOLS/')

from L_aln_tools import muscle_aln, trim_aln_gaps, cons_prof, add_consensus
from L_hist_aln import *
from L_seq_subset import *
from L_shade_aln import *
from L_plot4seq import *
from hist_ss import get_hist_ss_in_aln_for_html, get_hist_ss_in_aln_as_string
from L_aln2html import aln2html
from ete2 import NCBITaxa
from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node,AttrFace,TextFace
from Bio.Align import MultipleSeqAlignment
from math import log


ncbi = NCBITaxa()

def main():
    title=''
    #1. Getting data
    df=pd.read_csv('int_data/seqs_rs_redef.csv') #Histone types info
    fasta_dict=pickle.load( open( "int_data/fasta_dict.p", "rb" )) #Sequences
    # exit()
    
    #2. Filtering
    ##########
    #2.1. Narrow by variant/type
    title+='Canonical H2B'
    # f_df=df[(df['hist_var']=='canonical_H4')]
    # f_df['hist_var']='canonical_H4'
    f_df=df[((df['hist_var']=='canonical_H2B')|(df['hist_var']=='H2B.1'))]
    # &(df['taxid']!=118797)&(df['taxid']!=310752)
    # f_df=df[(df['hist_type']=='H2A')]
    # exit()
    print len(f_df)
    #2.2. #####select one variant per taxid
    # title+=' 1ptax'
    f_df=f_df.sort(['RefSeq'],ascending=False) # so that RefSeq record get priority on removing duplicates
    f_df=f_df.drop_duplicates(['taxid','hist_var'])


    # exit()
    #2.3. Filter by list of taxonomy clades   
    ################
    title+=' across cellular organisms'
    # parent_nodes=[9443] #131567 - cellular organisms, 7215 4930 Drosophila and yeast, 9443 - primates
    parent_nodes=[131567] #131567 - cellular organisms, 7215 4930 Drosophila and yeast, 9443 - primates
    #33682 - euglenozoa
    #6656 - arthropods
    # 4751 - fungi
    print "Selecting taxonomic subset"
    taxids=list(parent_nodes)
    for i in parent_nodes:
        taxids.extend(ncbi.get_descendant_taxa(i,intermediate_nodes=True))
    f_df=f_df[f_df['taxid'].isin(taxids)]
    print len(f_df)
    # exit()
    
    #2.4 Take one representative per specific taxonomic rank.
    ################
    title+=', one sequence per order'
    print "Pruning taxonomy"
    #Common ranks: superorder-order-suborder-infraorder-parvorder-superfamily-family-subfamily-genus-species-subspecies
    seqtaxids=list(f_df['taxid']) #old list
    new_seqtaxids=subsample_taxids(seqtaxids,rank='order') #new subsampled list
    f_df=f_df[f_df['taxid'].isin(new_seqtaxids)] #remake the dataframe
    # print "---"
    print len(f_df)
    # exit()


    #2.5. Check seq for sanity
    ################
    # title+=' seqQC '

    print "Checkig sequence quality"
    newgis=list()
    for i,row in f_df.iterrows():
        gi=row['gi']
        seq=fasta_dict[str(gi)].seq
        hist_type=row['hist_type']
        hist_var=row['hist_var']
        if(check_hist_length(seq,hist_type,hist_var,5)&check_hist_core_length(seq,hist_type,5)):
            newgis.append(gi)
    f_df=f_df[f_df['gi'].isin(newgis)] #remake the dataframe
    print len(f_df)
    # print list(f_df['gi'])
    # exit()

    #3. Make a list of seq with good ids and descriptions
    ####################
    f_fasta_dict={key: value for (key,value) in fasta_dict.iteritems() if int(key) in list(f_df['gi'])}
    print len(f_fasta_dict)
    taxid2name = ncbi.get_taxid_translator(list(f_df['taxid']))
    #Relabel sequences gi=> type and organism
    f_fasta_dict={key: SeqRecord(id=key, description=f_df.loc[f_df.gi==int(key),'hist_var'].values[0]+' '+taxid2name[f_df.loc[f_df.gi==int(key),'taxid'].values[0]],seq=value.seq) for (key,value) in f_fasta_dict.iteritems() }
    #with arbitrary index
    # f_fasta_dict_rel={key: SeqRecord(id=str(index), description=f_hist_df.loc[f_hist_df.gi==key,'hist_var'].values[0]+' '+taxid2names[f_hist_df.loc[f_hist_df.gi==key,'taxid'].values[0]],seq=f_fasta_dict[key].seq) for (index,key) in enumerate(f_fasta_dict) }
    # exit()

    #4. Make MSA
    #################
    #Here we construct MSA
    msa=muscle_aln(f_fasta_dict.values())
    AlignIO.write(msa, "results/h2b_ca_cellular.fasta", "fasta")

    msa_annot=MultipleSeqAlignment([SeqRecord(Seq(''.join(get_hist_ss_in_aln_as_string(msa)).replace(' ','-')),id='annotation',name='')])
    msa_annot.extend(msa)
    AlignIO.write(msa_annot, "results/h2b_ca_cellular_annot.fasta", "fasta")

    for i in range(len(msa)):
        gi=msa[i].id
        msa[i].description=f_fasta_dict[gi].description.replace('canonical','ca')
    msa.sort(key=lambda x: x.description)


    #5. Visualize MSA
    aln2html(msa,'results/h2b_ca_cellular.html',features=get_hist_ss_in_aln_for_html(msa,'H2B',0),title="canonical H2B in cellular organisms",description=True,field1w=10,field2w=35)


    #6. Trim alignment - this is optional
    #6.1. Trim gaps
    title+=', N-tail removed'
    # msa_tr=trim_aln_gaps(msa,threshold=0.8)

    #6.2. Trim to histone core sequence
    # msa_tr=trim_hist_aln_to_core(msa)
    msa_tr=msa[:,120:]

    #7. Vizualize MSA with ete2.
    taxid2gi={f_df.loc[f_df.gi==int(gi),'taxid'].values[0]:gi for gi in list(f_df['gi'])}
    gi2variant={gi:f_df.loc[f_df.gi==int(gi),'hist_var'].values[0] for gi in list(f_df['gi'])}

    msa_dict={i.id:i.seq for i in msa_tr}
    print taxid2gi
    t = ncbi.get_topology(list(f_df['taxid']),intermediate_nodes=False)
    a=t.add_child(name='annotation')
    a.add_feature('sci_name','annotation')
    t.sort_descendants(attr='sci_name')
    ts = TreeStyle()
    def layout(node):
        # print node.rank
        # print node.sci_name
        if getattr(node, "rank", None):
            if(node.rank in ['order','class','phylum','kingdom']):   
                rank_face = AttrFace("sci_name", fsize=7, fgcolor="indianred")
                node.add_face(rank_face, column=0, position="branch-top")
        if node.is_leaf():
            sciname_face = AttrFace("sci_name", fsize=9, fgcolor="steelblue")
            node.add_face(sciname_face, column=0, position="branch-right")
        if node.is_leaf() and not node.name=='annotation':
            s=str(msa_dict[str(taxid2gi[int(node.name)])])
            seqFace = SeqMotifFace(s,[[0,len(s), "seq", 10, 10, None, None, None]],scale_factor=1)
            add_face_to_node(seqFace, node, 0, position="aligned")
            gi=taxid2gi[int(node.name)]
            add_face_to_node(TextFace(' '+str(gi)+' '),node,column=1, position = "aligned")
            add_face_to_node(TextFace('      '+str(int(node.name))+' '),node,column=2, position = "aligned")
            add_face_to_node(TextFace('      '+str(gi2variant[gi])+' '),node,column=3, position = "aligned")

        if node.is_leaf() and node.name=='annotation':
            s=get_hist_ss_in_aln_as_string(msa_tr)
            seqFace = SeqMotifFace(s,[[0,len(s), "seq", 10, 10, None, None, None]],scale_factor=1)
            add_face_to_node(seqFace, node, 0, position="aligned")
            add_face_to_node(TextFace(' '+'NCBI_GI'+' '),node,column=1, position = "aligned")
            add_face_to_node(TextFace('       '+'NCBI_TAXID'+' '),node,column=2, position = "aligned")
            add_face_to_node(TextFace('       '+'Variant'+'       '),node,column=3, position = "aligned")



    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(title, fsize=20), column=0)
    t.render("results/h2b_ca_cellular.svg", w=6000, dpi=300, tree_style=ts)


    #10. Conservation
    features=get_hist_ss_in_aln_for_shade(msa_tr,below=True)
    cn=add_consensus(msa_tr,threshold=0.5)[-2:-1]

    # Below are three methods that we find useful.
    # plot_prof4seq('cons_sofp_psic',map(float,cons_prof(msa_tr,f=2,c=2)),cn,features,axis='conservation')
    plot_prof4seq('results/h2b_ca_cellular_cons_ent_unw',map(lambda x:log(20)+x,map(float,cons_prof(msa_tr,f=0,c=0))),cn,features,axis='conservation',title='Conservation, canonical H2B cellular organisms')
    
    # plot_prof4seq('cons_sofp_unw',map(float,cons_prof(msa_tr,f=0,c=2)),cn,features,axis='conservation')
    plot_prof4seq('results/h2b_ca_cellular_cons_sofp_unw_renorm1',map(float,cons_prof(msa_tr,f=0,c=2,m=1)),cn,features,axis='conservation',title='Conservation, canonical H2B cellular organisms')
    plot_prof4seq('results/h2b_ca_cellular_cons_sofp_psic_renorm1',map(float,cons_prof(msa_tr,f=2,c=2,m=1)),cn,features,axis='conservation',title='Conservation, canonical H2B cellular organisms')
    
    # plot_prof4seq('cons_ent_psic',map(lambda x:log(20)+x,map(float,cons_prof(msa_tr,f=2,c=0))),cn,features,axis='conservation')

   
#we get an alignment - we need to get conservation profile and visualize it on annotated histone sequence
    
    #11. Subfamily specific sites


    


    #12.Phylogenetic trees




if __name__ == '__main__':
    main()
    
            