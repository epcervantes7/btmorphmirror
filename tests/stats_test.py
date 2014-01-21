"""
Test routines for the btstructs.py file
"""

import sys

import numpy as np
import matplotlib.pyplot as plt

# import btstructs
# import btstats
import btmorph

swc_tree = None
stats = None

test_file_name = 'tests/v_e_moto1.CNG.swc'

def test_0_load_swc():
    """
    Test whether SWC files are correctly loaded
    """
    global swc_tree
    swc_tree = btmorph.STree2()
    swc_tree.read_SWC_tree_from_file(test_file_name)
    all_nodes = swc_tree.get_nodes()
    print 'len(all_nodes)', len(all_nodes)
    assert(len(all_nodes) == 562)

def test_0_load_stats() :
    """
    Check whether the statitics package can be loaded
    """
    global stats
    stats = btmorph.BTStats(swc_tree)

def test_global_bifurcations() :
    """
    Bifurcation count
    """
    no_bifurcations = stats.no_bifurcations()
    print 'no_bifurcations=%f' % (no_bifurcations)
    assert(no_bifurcations == 122)

def test_global_terminals() :    
    """
    Terminal point count
    """
    no_terminals = stats.no_terminals()
    print 'no_terminals=%f'  % (no_terminals)
    assert(no_terminals == 132)

def test_global_stems() :
    """
    Stem count
    """    
    no_stems = stats.no_stems()
    print 'no_stems=%s'  % (no_stems)
    assert(no_stems == 10)

def test_global_totallength() :    
    """
    Total length count
    """
    total_length = stats.total_length()
    print 'total length=%f' % (total_length)
    assert(78849 < total_length < 78850)

def test_global_somasurface() :    
    """
    Soma surface
    """
    soma_surface = stats.approx_soma()
    print 'soma surface=%f' % (soma_surface)
    assert(45238 < soma_surface < 45239)
    
def test_segment_length() :
    """
    Compute total length as sum of incoming segment lengths
    """
    bif_nodes = stats._bif_points
    term_nodes = stats._end_points
    all_nodes = bif_nodes +term_nodes
    total_length = 0
    for node in all_nodes :
        total_length = total_length + stats.get_segment_pathlength(node)
    print 'total_length=', total_length
    assert(78849 < total_length < 78850)

def test_terminal_lengths() :
    """
    Check terminal point lengths
    """
    term_path_lengths = []
    term_euclidean_lengths = []
    term_contractions = []
    for node in stats._end_points :
        term_path_lengths.append(stats.get_pathlength_to_root(node))
        term_euclidean_lengths.append(stats.get_Euclidean_length_to_root(node))
        term_contractions.append( term_euclidean_lengths[-1] / term_path_lengths[-1]  )
    print 'min/max path: %f - %f' % (min(term_path_lengths), max(term_path_lengths))
    print 'min/max euclid: %f - %f' % (min(term_euclidean_lengths), max(term_euclidean_lengths))
    print 'min/max contraction: %f - %f' % (min(term_contractions), max(term_contractions))
    assert(1531 < max(term_euclidean_lengths) <1532)
    assert(1817 < max(term_path_lengths) < 1819)

    
def test_degree() :
    """
    Topological degree
    """
    max_degree = stats.degree_of_node(swc_tree.get_root())
    print 'max_degree = ', max_degree
    assert(max_degree == 134)
    
def test_order() :
    """
    Topological order
    """
    max_order= -1
    min_order = 100000
    for node in swc_tree.get_nodes() : # skip the root
        order = stats.order_of_node(node)
        if order > max_order :
            max_order = order
        if order < min_order :
            min_order = order
    print 'min_order=%f, max_order=%f' % (min_order,max_order)
    assert(max_order == 9)

def test_partition_asymmetry():
    """
    Parition asymmetry
    """
    pa = []
    for node in stats._bif_points :
        pa.append(stats.partition_asymmetry(node))
    avg_pa = np.mean(pa)
    max_pa = max(pa)
    min_pa = min(pa)
    print 'avg_pa=%f, min_pa=%f, max_pa=%f' % (avg_pa, min_pa, max_pa)    
    assert(0.43 < avg_pa < 0.45)

def test_surface() :
    """
    Total neurite surface
    """
    total_surf, all_surfs = stats.total_surface()
    print 'total_surf=%f' % (total_surf)
    assert(512417 < total_surf < 512419)    

def test_volume() :
    """
    Total neurite volume
    """
    total_vol, all_vols = stats.total_volume()
    print 'total_volume=%f' % (total_vol)
    assert(390412 < total_vol < 390414)    

def ttest_bifurcation_sibling_ratio_local() :
    ratios = []
    for node in stats._bif_points :
        ratio = stats.bifurcation_sibling_ratio(node,where='local')
        ratios.append(ratio)
    print 'mean(ratios_local)=', np.mean(ratios)
    assert(1.31 < np.mean(ratios) < 1.32)


def ttest_bifurcation_sibling_ratio_remote() :
    ratios = []
    for node in stats._bif_points :
        ratio = stats.bifurcation_sibling_ratio(node,where='remote')
        ratios.append(ratio)
    print 'mean(ratios_remote)=', np.mean(ratios)    
    assert(1.16 < np.mean(ratios) < 1.17)

# def test_bifurcation_amplitude_local() :
#     all_ampl = []
#     for node in stats._bif_points :
#         ampl = stats.bifurcation_angle(node, where='local')[0]
#         all_ampl.append(ampl)
#     print 'min=%f max(ample)=%f, mean(ampl)=%f' % (np.min(all_ampl),np.max(all_ampl),np.mean(all_ampl))
#     import matplotlib.pyplot as plt
#     plt.figure(1)
#     plt.hist(all_ampl,color='blue',alpha=0.5)

# def test_bifurcation_amplitude_remote() :
#     all_ampl = []
#     for node in stats._bif_points :
#         ampl = stats.bifurcation_angle(node, where='remote')[0]
#         all_ampl.append(ampl)
#     print 'min=%f max(ample)=%f, mean(ampl)=%f' % (np.min(all_ampl),np.max(all_ampl),np.mean(all_ampl))
#     import matplotlib.pyplot as plt
#     plt.figure(1)
#     plt.hist(all_ampl,color='red',alpha=0.5)
#     # plt.show()
    
def ttest_ralls_ratio() :
    """
    Binary search for rall's power
    """
    all_p = []
    for node in stats._bif_points :
        p = stats.bifurcation_ralls_ratio(node) 
        all_p.append(p)
        print node, '-> p=', p
    all_p = np.array(all_p)
    all_pp = []
    for n in all_p :
        if not np.isnan(n) :
            all_pp.append(n)
    print 'min_p=%f,avg_p=%f, max_p=%f' % (np.min(all_pp),np.mean(all_pp),np.max(all_pp))
    # p = stats.bifurcation_ralls_ratio(stats._bif_points[1])
    plt.hist(all_pp)
    plt.show()

def ttest_ralls_ratio2() :
    """
    Binary search for rall's power
    """
    all_p = []
    for node in stats._bif_points :
        p = stats.bifurcation_ralls_ratio2(node) 
        all_p.append(p)
        print node, '-> p=', p
    all_p = np.array(all_p)
    all_pp = []
    for n in all_p :
        if not np.isnan(n) :
            all_pp.append(n)
    print 'min_p=%f,avg_p=%f media=%f, max_p=%f' % (np.min(all_pp),np.mean(all_pp),np.median(all_pp),np.max(all_pp))
    # p = stats.bifurcation_ralls_ratio(stats._bif_points[1])
    plt.hist(all_pp)
    plt.show()
    
# def test_plot_SWC_3D() :
#     """
#     Plot SWC file in 3D, saves figures in tests/test_figure.pdf
#     """
#     import btviz
#     import matplotlib.pyplot as plt
#     btviz.plot_3D_SWC(file_name=test_file_name)
#     plt.savefig('tests/test_figure3D.pdf')
#     assert(True) # if not crashed, it's ok
    
# def test_plot_SWC_2D() :
#     """
#     Plot SWC file in 3D, saves figures in tests/test_figure.pdf
#     """
#     import btviz
#     import matplotlib.pyplot as plt
#     btviz.plot_2D_SWC(file_name=test_file_name)
#     plt.savefig('tests/test_figure2D.pdf')    
#     assert(True) # if not crashed, it's ok

# def test_dendrogram() :
#     import btviz
#     import matplotlib.pyplot as plt
#     btviz.plot_dendrogram(file_name=test_file_name)
#     plt.savefig('tests/test_figure_dendrogram.pdf')
#     assert(True) # if not crashed, it's ok
