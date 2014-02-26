import os
import sys
import glob
import cPickle as pickle
import matplotlib.pyplot as plt

import btmorph

def perform_1D_population_analysis(destination):
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
        # plt.clf()
        # btmorph.plot_2D_SWC(f,outN=cell_name+"_2D.pdf")
        # plt.clf()
        # btmorph.plot_3D_SWC(f,outN=cell_name+"_3D.pdf")    

    """ 1D features, without dependencies on other quantities
        Both global (1 per neuron) and local (N per neuron) features"""

    one_d_stats = {}

    # bulk statistics of global features
    funcs = ["total_length","no_bifurcations","no_terminals","no_stems",\
             "total_volume","total_surface"]
    for cell_name in individual_stats:
        for func in funcs:
            method = getattr(individual_stats[cell_name], func)
            ret = method()
            if isinstance(ret,tuple):
                ret = ret[0]
            if func in one_d_stats:
                one_d_stats[func].append(ret)
            else:
                one_d_stats[func] = []
                one_d_stats[func].append(ret)

    # extract the max order            
    for cell_name in individual_stats:
        max_order = 0
        for node in swc_trees[cell_name].get_nodes():
            order = individual_stats[cell_name].order_of_node(node)
            if order > max_order:
                max_order = order
        if "max_order" in one_d_stats:
            one_d_stats["max_order"].append(max_order)
        else:
            one_d_stats["max_order"] = []
            one_d_stats["max_order"].append(max_order)

    # extract the "inter-bifurcation segment length"        
    for cell_name in individual_stats:
        bif_nodes = individual_stats[cell_name]._bif_points
        for node in bif_nodes:
            L = individual_stats[cell_name].get_segment_pathlength(node)
            if "inter_bif_L" in one_d_stats:
                one_d_stats["inter_bif_L"].append(L)
            else:
                one_d_stats["inter_bif_L"] = []
                one_d_stats["inter_bif_L"].append(L)

    # extract the terminal segment lengths
    # extract the path length / euclidean distance to terminal tips
    for cell_name in individual_stats:
        term_nodes = individual_stats[cell_name]._end_points
        for node in term_nodes:
            L = individual_stats[cell_name].get_segment_pathlength(node)
            path_L = individual_stats[cell_name].get_pathlength_to_root(node)
            eucl_L = individual_stats[cell_name].get_Euclidean_length_to_root(node)
            if "term_L" in one_d_stats:
                one_d_stats["term_L"].append(L)
            else:
                one_d_stats["term_L"] = []
                one_d_stats["term_L"].append(L)

            if "term_path_L" in one_d_stats:
                one_d_stats["term_path_L"].append(path_L)
            else:
                one_d_stats["term_path_L"] = []
                one_d_stats["term_path_L"].append(path_L)

            if "term_eucl_L" in one_d_stats:
                one_d_stats["term_eucl_L"].append(eucl_L)
            else:
                one_d_stats["term_eucl_L"] = []
                one_d_stats["term_eucl_L"].append(eucl_L)


    for func in one_d_stats:
        plt.figure(0)
        plt.clf()
        plt.hist(one_d_stats[func])
        plt.xlabel(func)
        plt.ylabel("#")
        plt.savefig("pop_1D_"+func+"_hist.pdf")
        p_name = "pop_1D_"+func+"_raw.pkl"
        pickle.dump(one_d_stats[func], open(p_name,"w"))

if __name__=="__main__":
    destination = "."
    if len(sys.argv) == 2:
        destination = sys.argv[1]
    perform_1D_population_analysis(destination)

    
