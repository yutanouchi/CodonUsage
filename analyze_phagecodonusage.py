#!/usr/bin/python


from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt
import matplotlib
from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist
import codonanalyzer


filepath='sequence_data/lambda.gb'
record=SeqIO.parse(open(filepath,'rU'),'genbank').next()
book=codonanalyzer.CodonBook()
book.add_genes(record)

# import ipdb; ipdb.set_trace()#

df=book.normalize_by_totalaa().ix[:,1:].transpose()
df=df.ix[:,book.regulatorycodons]

# Compute and plot first dendrogram.
fig = plt.figure(figsize=(18,12))
ax1 = fig.add_axes([0.05,0.1,0.2,0.8])
y=linkage(df,method='ward')
z=dendrogram(y,orientation='right')
ax1.set_xticks([])
ax1.set_yticks([])


# Plot distance matrix.
axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
idx1 = z['leaves']
# idx2 = Z2['leaves']
df = df.ix[idx1,:]
# D = D.ix[:,idx2]
im = axmatrix.matshow(df, aspect='auto', origin='lower', cmap=matplotlib.cm.YlOrRd)
axmatrix.set_xticks(range(len(df.columns)))
axmatrix.set_yticks(range(len(df.index)))
axmatrix.set_xticklabels(df.columns,rotation=90)
axmatrix.set_yticklabels(df.index.map(lambda x: book.lookup_locustag(x)))

# Plot colorbar.
axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
plt.colorbar(im,cax=axcolor)
fig.savefig('cluster.pdf')
fig.show()




