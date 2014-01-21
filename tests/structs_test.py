'''
Test routines for the btstructs2.py file
'''

import sys
#sys.path.append('..')

#from btstructs import STree, SNode,P3D
#import btstructs2

import btmorph
from btmorph import STree2

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

