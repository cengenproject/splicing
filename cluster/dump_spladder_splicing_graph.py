import pickle
import os
import numpy as np        

# Note: these commands were executed interactively

gene_file='/path_to/spladder/spladder_AB_bss2/spladder/genes_graph_conf3.merge_graphs.validated.pickle'

(genes, events) = pickle.load(open(gene_file, 'rb'), encoding='latin1')

len(genes) # save for later, should be ~46,910

# if necessary to avoid appending
# os.remove('vertices_AVE.csv')

f=open('200726_vertices.csv','a')

for gene in genes:
	f.write('> '+gene.name+' '+gene.chr+' '+gene.strand+'\n')
	np.savetxt(f,gene.splicegraph.vertices, fmt="%u, "*gene.splicegraph.vertices.shape[1])

f.close()




# Check for mbk-2, should have that exon at 13012366-13012566
for gene in genes:
    if gene.name == "WBGene00003150":
        print("Found it")
        my_gene = gene