'''
Test routines for the btstructs2.py file
'''

import sys
import numpy
import random as rndm
import itertools
#sys.path.append('..')

#from btstructs import STree, SNode,P3D
#import btstructs2

import btmorph
from btmorph import STree2
from btmorph import VoxelGrid

def test_load_swc():
    '''
    Test whether SWC files are correctly loaded
    '''
    swc_tree = btmorph.STree2()
    swc_tree.read_SWC_tree_from_file('tests/v_e_moto1.CNG.swc')
    all_nodes1 = swc_tree.get_nodes()

    print '\nlen(swc_tree1)', len(all_nodes1) 

    assert(len(all_nodes1) == 562)

def test_load_and_write_swc():
    '''
    Test whether SWC treesare correctly written to file
    '''
    swc_tree = STree2()
    swc_tree.read_SWC_tree_from_file('tests/v_e_moto1.CNG.swc')
    swc_tree.write_SWC_tree_to_file('tests/moto_1_outputted.swc')
    swc_tree2 = STree2()
    swc_tree2.read_SWC_tree_from_file('tests/moto_1_outputted.swc')
    print 'len(swc_tree2)', len(swc_tree2.get_nodes()) 
    
    assert(len(swc_tree2.get_nodes()) == 562)
    

if __name__ == '__main__' :
    test_load_swc()

def test_voxelGrid_adjust_dimensions():
    """
    Test whether x:y:z = rx:ry:rz, where (x,y,z) - output dimensions, (rx,ry,rz) - resolution.
        And none of the original dimensions can be bigger than the output (expansion only)
    Special cases: One dimension and corresponding resolution are zero => no modification needed
                   More than one dimension is zero => return nothing
                   One or more dimension/resolution is zero while corresponding resolution/dimension is not => return nothing
    """   
    epsilon = 0.001
    # One dimension and corresponding resolution are zero
    for i in range(0, 3):
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)] #[rx, ry, rz]        
        #One of the dimensions is zero
        dimensions[i] = 0
        resolution[i] = 0
        # No need to expand anything
        #[x,y,z] = dimensions
        #[rx,ry,rz] = VoxelGrid.adjustDimensions(dimensions, resolution)
        #assert(abs(x - rx) < epsilon)
        #assert(abs(y - ry) < epsilon)
        #assert(abs(z - rz) < epsilon)
    # More than one dimension is zero => return nothing
    dimensions = numpy.array([0, 0, 0])
    resolution = [0, 0, 0]
    assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    for i in range(0, 3):
        dimensions = numpy.array([0, 0, 0])
        resolution = [0, 0, 0]
        #One of the dimensions is not zero
        dimensions[i] = rndm.randint(1000, 2000)/10.0
        resolution[i] = rndm.randint(1000, 2000)
        assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    # One or more dimension/resolution is zero while corresponding resolution/dimension is not => return nothing
    # All
    dimensions = numpy.array([0, 0, 0])
    resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]
    assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
    resolution = [0, 0, 0]
    assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    # 1 (2 are taken care of in earlier)
    for i in range(0,3):
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)] #[rx, ry, rz]        
        dimensions[i] = 0
        assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)] #[rx, ry, rz]        
        resolution[i] = 0
        assert(VoxelGrid.adjustDimensions(dimensions, resolution) == None)
    # Everything is not zero
    # Random floating point from 100 to 200 with 1 digit after comma
    for i in range(0,10):
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)] #[rx, ry, rz]
        print(dimensions)
        print(resolution)    
        [x,y,z] = dimensions
        [rx,ry,rz]=resolution
        [nx,ny,nz] = VoxelGrid.adjustDimensions(dimensions, resolution)
        # None of the original dimensions can be bigger than the output (expansion only)
        assert(x <= nx) #x
        assert(y <= ny) #y
        assert(z <= nz) #z
        # x:y:z =?= rx:ry:rz
        print(abs(nx/ny - float(rx)/float(ry)))
        print(abs(ny/nz - float(ry)/float(rz)))
        assert(abs(nx/ny - float(rx)/float(ry)) < epsilon)
        assert(abs(ny/nz - float(ry)/float(rz)) < epsilon)    
    
    
    