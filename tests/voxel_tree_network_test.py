# -*- coding: utf-8 -*-
"""
Created on Fri Jul 11 18:21:38 2014
Test routines for VoxelizedTreeNetwork class
@author: Irina Reshodko
"""
import btmorph
from btmorph import VoxelGrid
from btmorph import plot_3D_SWC
from btmorph import pca_project_tree
from btmorph import VoxelizedTreeNetwork
from  nose.tools import raises

def prepare_trees():
    tree1 = btmorph.STree2().read_SWC_tree_from_file("tests/CTh5080306F.CNG.swc")#horton-strahler_test_wiki_3pointsoma.swc")
    tree1 = pca_project_tree(tree1)
    tree2 = btmorph.STree2().read_SWC_tree_from_file("tests/moto_1_outputted.swc")
    tree2 = pca_project_tree(tree2) 
    return tree1, tree2
    
def test_init_ok():
    """
    Test if VoxelizedTreeNetwork initializes properly: normal case
    """
    tree1, tree2 =  prepare_trees()
    # 3D
    dim = (1500, 1500, 260)
    voxelSz = 10.0
    percInd = 1
    network = VoxelizedTreeNetwork(dim, voxelSz, [tree1, tree2], percInd)
    network.vg.plot()
    assert(len(network.vg.grid) > 0)
    res = network.vg.res
    t = res[percInd]
    res[percInd] = res[0]
    res[0] = t
    assert(res == network.perc.size)
    # 2D
    dim = (1500, 1500, 0)
    network = VoxelizedTreeNetwork(dim, voxelSz, [tree1, tree2])
    network.vg.plot()
    print btmorph.BTStats(tree2).total_dimension()
    print btmorph.BTStats(tree1).total_dimension()
    print network.vg
    assert(len(network.vg.grid) > 0)

@raises(TypeError)
def init_bad_argument(dim, vsz, trees):
    """
    Test if TypeError is thrown in init if an argument is invalid
    """
    network = VoxelizedTreeNetwork(dim, vsz, trees)

def test_init_bad_arguments():
    """
    Test if init reacts properly to invalid argument
    """
    bad_dims = [[], [1,2], [1,3,5,6], [1], [-100, 100, 100], [100, -100, 100], [100, 100, -100]]
    bad_vsz = [0, -100]
    for dim in bad_dims:
        yield init_bad_argument, dim, 1, None
    for vsz in bad_vsz:
        yield init_bad_argument, [100, 100, 100], vsz, None
        
def test_add_tree():
    """
    Test if a tree is added properly to the network
    """
    dim = (1500, 1500, 260)
    voxelSz = 10.0
    tree1, tree2 =  prepare_trees()
    network = VoxelizedTreeNetwork(dim, voxelSz)
    network.add_tree(tree1)
    sz1 = len(network.vg.grid)
    assert(sz1 > 0)
    network.add_tree(tree2)
    sz2 = len(network.vg.grid)
    assert(sz2 > sz1)
    network.add_tree(None)
    assert(len(network.vg.grid) == sz2)
    
    

    