#!/Users/epcervantes/anaconda/bin/python
import numpy
import matplotlib.pyplot as plt

import sys

sys.path.insert(0, '../btmorph')
import btmorph
swc_tree= btmorph.STree2()
types=range(1,10)
swc_tree.read_SWC_tree_from_file("examples/data/v_e_moto1.CNG.swc",types)

stats = btmorph.BTStats(swc_tree)

# get the total length
total_length = stats.total_length()
print "total_length = %f" % total_length

# get the max degree, i.e., degree of the soma
max_degree = stats.degree_of_node(swc_tree.get_root())

# generate and save the dendrogram
btmorph.plot_dendrogram("examples/data/v_e_moto1.CNG.swc")
plt.savefig('examplar_dendrogram.pdf')