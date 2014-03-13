import os
import sys
import glob
import cPickle as pickle
import matplotlib.pyplot as plt

import btmorph

def filter_and_save_SWC(destination,types=range(10),prefix="_filtered"):
    """
    Removes points from a SWC structure and saves the new SWC to a file.

    Can be used to remove unwanted structures that are identifiable by \
    the type-field in the SWC description.

    Parameters
    -----------
    destination : string
        string with the location of where to find the SWC files.
    types : list of int
        types that are to be remained in the SWC file.
    """
    # change to a directory of choice for I/O
    os.chdir(destination)    

    # load morphologies and initialize statistics
    all_f = glob.glob("*.swc")

    for f in all_f:
        temp_tree = btmorph.STree2()
        temp_tree.read_SWC_tree_from_file(f,types=types)
        outN = f.split(".swc")[0]+prefix+".swc"
        temp_tree.write_SWC_tree_to_file(outN)
