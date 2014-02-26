import os
import sys
import glob
import cPickle as pickle
import matplotlib.pyplot as plt

import btmorph

def _get_node_features(stats,node,term=False):
    O = stats.order_of_node(node)
    D = stats.order_of_node(node)
    ED = stats.get_Euclidean_length_to_root(node)
    PL = stats.get_pathlength_to_root(node)
    SPL = stats.get_segment_pathlength(node)
    SED = stats.get_segment_Euclidean_length(node)
    if term:
        return O,D,ED,PL,SPL,SED
    else:
        PA = stats.partition_asymmetry(node)
        return O,D,ED,PL,SPL,SED,PA

def perform_2D_analysis(destination):
    """function description"""

    # change to a directory of choice for I/O
    os.chdir(destination)
    
    # load morphologies and initialize statistics
    all_f = glob.glob("*.swc")
    swc_trees = {}
    individual_stats = {}
    for f in all_f:
        print "f: ", f
        cell_name = f.split(".")[0]
        temp_tree = btmorph.STree2()
        temp_tree.read_SWC_tree_from_file(f)
        swc_trees[cell_name] = temp_tree
        temp_stats = btmorph.BTStats(temp_tree)
        individual_stats[cell_name] = temp_stats

    """ 2D features calculated per neuron and for the whole population.
        For each bifurcation and terminal point, the following features \
        at that location will be recorded: order, degree, asymmetry_index, \
        path length, euclidean distance,(bifurcation angle, local)"""

    for cell_name in individual_stats:
        print "analyzing cell %s" % cell_name
        bif_nodes = individual_stats[cell_name]._bif_points
        term_nodes = individual_stats[cell_name]._end_points

        bif_to_write = ""
        for node in bif_nodes:
            O,D,ED,PL,SPL,SED,PA= \
              _get_node_features(individual_stats[cell_name],node)
            bif_to_write += "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n" % (O,D,ED,PL,SPL,SED,PA)
        bif_writer = open(cell_name+"_bifs_2D.txt","w")
        bif_writer.write(bif_to_write)
        bif_writer.flush()
        bif_writer.close()

        term_to_write = ""
        for node in term_nodes:
            O,D,ED,PL,SPL,SED= \
              _get_node_features(individual_stats[cell_name],node,term=True)
            term_to_write += "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n" % (O,D,ED,PL,SPL,SED)
            
        term_writer = open(cell_name+"_terms_2D.txt","w")
        term_writer.write(term_to_write)
        term_writer.flush()
        term_writer.close()

if __name__=="__main__":
    destination = "."
    if len(sys.argv) == 2:
        destination = sys.argv[1]
    perform_2D_analysis(destination)
