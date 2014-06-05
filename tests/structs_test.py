'''
Test routines for the btstructs2.py file
'''

import sys
import numpy
import random as rndm
from  nose.tools import raises
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

@raises(TypeError)
def VoxelGrid_adjustDim_junk(d, r):
    """
    Test if adjustDimensions raises TypeError if fed with junk
    """
    VoxelGrid.adjustDimensions(d, r)

@raises(IndexError, Exception)
def VoxelGrid_adjustDim_negative(d, r):
    """
    Test if adjustDimensions raises IndexError if fed with negative value
    """
    VoxelGrid.adjustDimensions(d, r)

def test_voxelGrid_adjust_dimensions():
    """
    Test whether x:y:z = rx:ry:rz, where (x,y,z) - output dimensions, (rx,ry,rz) - resolution.
        And none of the original dimensions can be bigger than the output (expansion only)
    Special cases: One dimension and corresponding resolution are zero => no modification needed
                   More than one dimension is zero => return nothing
                   One or more dimension/resolution is zero while corresponding resolution/dimension is not => return nothing
    """   
    epsilon = 0.001
    # Should not accept junk and values less than zero
    junk = [None, [], [4, 't', 4], "str"]
    dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
    resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]    
    for el in junk:
        yield VoxelGrid_adjustDim_junk, dimensions, junk
        yield VoxelGrid_adjustDim_junk, junk, resolution
    for i in range(0,3):
        dimensions = numpy.array([rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0, rndm.randint(1000, 2000)/10.0]) #[x, y, z]
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]    
        resolution[i] = -100
        yield VoxelGrid_adjustDim_negative, dimensions, resolution
        resolution = [rndm.randint(1000, 2000), rndm.randint(1000, 2000), rndm.randint(1000, 2000)]    
        dimensions[i] = -100
        yield VoxelGrid_adjustDim_negative, dimensions, resolution
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
    for i in range(0,10):
        # Random floating point from 100 to 200 with 1 digit after comma
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
    
def test_isPowerOfTwo() :
    """
    Tests VoxelGrid.isPowerOfTwo function
    Must return True for 2^m and False otherwise
    """
    assert(VoxelGrid.isPowerOfTwo(-100) == False)
    assert(VoxelGrid.isPowerOfTwo(120.56) == False)
    for i in range(0,32):
        res_true = VoxelGrid.isPowerOfTwo(pow(2, i))
        res_false = VoxelGrid.isPowerOfTwo(pow(2, i)*3)
        assert(res_true == True)
        assert(res_false == False)
        
@raises(TypeError)     
def voxelGrid_check_key_bad_type(key):
    """
    The key must be a tuple of tree integers greater or equal than zero and 
    less than the cossesponding resolution
    Here testing: bad type
    """
    VoxelGrid.checkKey([100,100, 100], key)
    
@raises(IndexError)
def voxelGrid_check_key_out_of_boundary(key):
    """
    The key must be a tuple of tree integers greater or equal than zero and 
    less than the cossesponding resolution
    Here testing: out of boundary
    """
    VoxelGrid.checkKey([200,200,200], key)

def test_gen_voxelgrid_checkKey():
    """
    The key must be a tuple of tree integers greater or equal than zero and 
    less than the cossesponding resolution
    """    
    bad_type = [100.4, ['r', 'r'], "str", None, (10.8, 30, 190), (40, 15.7, 90), (78, 15, 87.9)]
    bad_num = [(100, -1, 100), (-1, 100, 100), (100, 100, -1)]
    bad_index = [(100, 1000, 100), (1000, 100, 100), (100, 100, 1000)]
    good = [(100, 120, 40), (11, 33, 90)]
    #bad type 
    for el in bad_type:
        yield voxelGrid_check_key_bad_type, el
    #bad number
    for el in bad_num:        
        yield voxelGrid_check_key_out_of_boundary, el
    #out of boundary
    for el in bad_index:
        yield voxelGrid_check_key_out_of_boundary, el
    #good
    for el in good:
        assert(VoxelGrid.checkKey([200,200,200], el) == True)

@raises(IndexError)
def VoxelGrid_init_not_pow2(dims, res):
    """
    Test if VoxelGrid initialization raises an exception if not power of two
    """
    VoxelGrid(dims, res)

@raises(TypeError)
def VoxelGrid_init_wrongType(dims, res):
    """
    Test if VoxelGrid initialization raises an exception if fed with junk
    """
    VoxelGrid(dims, res)

def test_VoxelGrid_init():
    """
    Test if VoxelGrid initialization was successfull
    """
    epsilon = 0.001
    #Should not accept junk
    junk = [None, [], [4, 't', 4], "str"]
    for el in junk:
        yield VoxelGrid_init_wrongType, junk, junk
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    # Should accept only power of two as resolution
    for i in range(0,3):
        resolution = [1024, 1024, 1024] #[rx, ry, rz]
        resolution[i] = 5*1024
        print (dimensions)
        print (resolution)
        print('----')
        yield VoxelGrid_init_not_pow2, dimensions, resolution        
    resolution = [64, 128, 64] #[rx, ry, rz]
    vg = VoxelGrid(dimensions, resolution)
    good_dim = VoxelGrid.adjustDimensions(dimensions, resolution)
    #check dimensions
    for i in range(0, 3):
        assert(abs(good_dim[i] - vg.dim[i]) < epsilon)
    #check grid - should be all blank
    for i in range(0, resolution[0]):
        for j in range(0, resolution[1]):
            for k in range(0, resolution[2]):
                assert(vg[(i,j,k)] == False)

@raises(TypeError)
def VoxelGrid_setitem_badtype(vg, tpl):
    """
    VoxelGrid.setitem method should raise an error if the type of the key is invalid
    """
    vg[tpl] = True

@raises(TypeError)
def VoxelGrid_setitem_badassign(vg, tpl, val):
    """
    VoxelGrid.setitem method should raise an error if the type of the value to be assigned is not boolean
    """
    vg[tpl] = val

@raises(IndexError)
def VoxelGrid_setitem_badindex(vg, tpl):
    """
    VoxelGrid.setitem method should raise an error if index is out of range
    """
    vg[tpl] = True
    
@raises(TypeError)
def VoxelGrid_getitem_badtype(vg, tpl):
    """
    VoxelGrid.getitem method should raise an error if the type of the key is invalid
    """
    vg[tpl]

@raises(IndexError)
def VoxelGrid_getitem_badindex(vg, tpl):
    """
    VoxelGrid.getitem method should raise an error if index is out of range
    """
    vg[tpl]

def test_VoxelGrid_setitem():
    """
    Test if VoxelGrid setitem method works properly
    """
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    resolution = [64, 128, 64]
    vg = VoxelGrid(dimensions, resolution)    
    # Bad type
    bad_type = [1, 20.5, (1,2), ["s"], ("1", 2, 3), (1, "2", 3), (1, 2, "3"), "str", None]
    for el in bad_type:
        yield VoxelGrid_setitem_badtype, vg, el
    # bad index
    for i in range(0, 3):
        key = [30, 100, 60]
        key[i] = -40
        yield VoxelGrid_setitem_badindex, vg, tuple(key)
    key = (30, 100, 60)
    # Value to be assigned is not boolean
    for el in bad_type:
        yield VoxelGrid_setitem_badassign, vg, key, el
     # Good example
    vg[key] = True            
    assert(key in vg.grid and vg.grid[key] == True)
    for i in range(0, resolution[0]):
        for j in range(0, resolution[1]):
            for k in range(0, resolution[2]):
                if i == key[0] and j == key[1] and k == key[2]:
                    continue
                assert(not (i,j,k) in vg.grid)
    vg[key] = False
    assert(len(vg.grid) == 0)
    
def test_VoxelGrid_getitem():
    """
    Test if VoxelGrid getitem method works properly
    """
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    resolution = [64, 128, 64]
    vg = VoxelGrid(dimensions, resolution)    
    # Bad type
    bad_type = [1, 20.5, (1,2), ["s"], ("1", 2, 3), (1, "2", 3), (1, 2, "3"), "str", None]
    for el in bad_type:
        yield VoxelGrid_getitem_badtype, vg, el
    # bad index
    for i in range(0, 3):
        key = [30, 100, 60]
        key[i] = -40
        yield VoxelGrid_getitem_badindex, vg, tuple(key)
    # Good example
    key = (30, 100, 60)
    vg.grid[key] = True
    r = vg[key]
    assert(r == True)
    for i in range(0, resolution[0]):
        for j in range(0, resolution[1]):
            for k in range(0, resolution[2]):
                if i == key[0] and j == key[1] and k == key[2]:
                    continue
                print (key)
                assert(vg[(i,j,k)] == False)

def test_VoxelGrid_getsetitem():
    """
    Test if setitem/getitem combination is correct for VoxelGrid
    """
    # Random floating point from 100 to 200 with 1 digit after comma
    dimensions = numpy.array([rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0, rndm.randint(100, 200)/10.0]) #[x, y, z]
    resolution = [64, 128, 64]
    vg = VoxelGrid(dimensions, resolution)    
    key = (30, 100, 60)
    vg[key] = True
    assert(vg[key] == True)
    