"""
Created on Fri Jul 11 17:32:44 2014
@author: Irina Reshodko

Defines data structure to store a network of voxelized trees
"""
import math
from btmorph import VoxelGrid
from btmorph import Percolation
import sets

class VoxelizedTreeNetwork:
    """
    Defines data structure to store a network of voxelized trees
    """
    def __init__(self, dim, voxelSz, trees = None, percDirection = 0):
        """
        Initialize
        
        Parameters:
        -------
        dim : iterable
        Global dimensions of the network in x, y and z direction in micrometers.
        |  If dim[z] == 0 then a 2d voxel grid will be generated
        voxelSz : number
        Dimensions of one voxel in x, y and z direction in micrometers. Must be greater than 0 
        |  This parameter determines resolution of the voxel grid.
        |  Will be adjusted according to the requirements of the grid size (power of two)
        trees : iterable of :class:`btmorph.btstructs2.STree2`
        Trees to be added to the network. Can be added later one be one
        percDirection : int
        Integer >= 0 and <=2. Specifies direction of percolation. Can only be set once.
        """
        if not len(dim) == 3:
            raise TypeError('Dimensions: Iterable of length 3 expected, got ' + type(dim) + "of length " + len(dim))
        if voxelSz <= 0:
            raise TypeError('Voxel size must be greater than zero. Got ' + str(voxelSz))
        for el in dim:
            if el < 0:
                raise TypeError('Dimensions must be greater than or equal to zero. Got ' + str(dim))
        self.dim = dim
        dx,dy,dz = self.dim
        res = [int(2**round(math.log(dx/voxelSz, 2))), int(2**round(math.log(dy/voxelSz, 2))), 1]
        if dz > 0.0:
            res[2] = int(2**round(math.log(dz/voxelSz, 2)))
        self.trees = []
        self.vg = VoxelGrid(self.dim, res)
        self.__percDirection__ = percDirection
        # Put percolation index in first place
        percRes = list(self.vg.res)
        t = percRes[self.__percDirection__]
        percRes[self.__percDirection__] = percRes[0]
        percRes[0] = t
        self.perc = Percolation(percRes)
        if trees == None:
            return
        for tree in trees:
            self.add_tree(tree)
        
    def add_tree(self, tree, offset = None):
        """
        Add a tree to the network
        
        Parameters:
        -------
        tree : :class:`btmorph.btstructs2.STree2`
        A tree to be voxelized and added to the network
        offset : iterable
        Tree position in the grid
        """
        if tree == None:
            return
        self.trees.append(tree)
        prevKeys = set(self.vg.grid.keys())
        self.vg.add_tree(tree, offset)
        newKeys = set(self.vg.grid.keys())
        diff = newKeys - prevKeys
        print self.vg
        for key in diff:
            # Put percolation index in first place
            key = map(lambda x: x+1, list(key))
            t = key[self.__percDirection__]
            key[self.__percDirection__] = key[0]
            key[0] = t
            self.perc.open(key)
    
