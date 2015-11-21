#!/usr/bin/python


from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist
import numpy as np
import codonanalyzer

plt.style.use('ggplot')

# load data
colifilepath='sequence_data/MG1655.gb'
colirecord=SeqIO.parse(open(colifilepath,'rU'),'genbank').next()
colibook=codonanalyzer.CodonBook()
colibook.add_genes(colirecord)
lambdafilepath='sequence_data/lambda.gb'
lambdarecord=SeqIO.parse(open(lambdafilepath,'rU'),'genbank').next()

# load expression data
proteincounts=pd.read_csv('lit_data/ProteinCounts_Li_2014.csv',index_col='Gene')
demand=colibook.codon_expression(proteincounts['MOPS complete'])
demand=demand/demand.max()

# calculate tRNA supply
coli_te=codonanalyzer.TransEff()
trna_abundance=pd.read_csv('lit_data/tRNA_Dong.csv',index_col=0)
supply=coli_te.weights_from_tRNAabundance(trna_abundance)

# translation efficiency index nte
nte=supply/demand
nte=nte/nte.max()

# phage codon usage
book=codonanalyzer.CodonBook()
book.add_genes(lambdarecord)
codonusage=book.normalize_by_totalaa().ix[:,1:]
gene_ai=np.exp(codonusage.apply(lambda x: x*np.log(nte)).sum())
import ipdb; ipdb.set_trace()#



# plotting
# plot_df=codon_tf.transpose().fillna(0)
plot_ds=gene_ai.sort_values().fillna(0)
fig, ax=plt.subplots(figsize=(18,12))
plot_ds.plot(kind='bar')
ax.set_xticklabels(plot_ds.index.map(lambda x:book.lookup_locustag(x)))

# fig.savefig('cluster.pdf')
fig.show()





import ipdb; ipdb.set_trace()#




