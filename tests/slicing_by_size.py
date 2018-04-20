#!/Users/epcervantes/anaconda/bin/python
##!/usr/bin/python
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore");
    import matplotlib.pyplot as plt
import numpy as np

import sys
import sys

sys.path.insert(0, '../btmorph')
import btmorph
import os

dictionary = {
    #'directory_name': class id
    'test': 0,
    'test1': 1,
    # 'chimpanzee_principal_pyramidal': 1,
    # 'human_principal_pyramidal': 2,
    # 'mouse_principal_ganglion': 3,
    # 'mouse_principal_pyramidal': 4,
    # 'rat_interneuron_gabaergic': 5,
    # 'rat_interneuron_nitrergic':6,
    # 'rat_principal_pyramidal_hippocampus': 7,
    # 'rat_principal_pyramidal_neocortex': 8,
    # 'various_principal_granule': 9,
    # 'various_principal_medium_spiny':10
}


if len(sys.argv) == 4:
    directory = "tests/" + str(sys.argv[1])
    slice_size = int(str(sys.argv[3]))
    new_directory = "tests/" + str(sys.argv[2]+"_" + str(slice_size))

    class_number=(dictionary[sys.argv[1]])
    from_level = 0
    to_level = 5
    counter = 0
    for i in range(from_level,to_level): # runs until the final level
        slice_begin = i
        slice_end = i + slice_size
        new_file_name = str(slice_begin)

        print(str(slice_begin) + '-' + str(slice_end))
        for fileName in os.listdir(directory): #for each file in the directory
            if fileName.endswith(".swc"):
                counter += 1
                swc_tree = btmorph.STree2()
                filePath = str(directory) + '/' + str(fileName)
                types = [1, 3, 4, 5, 6, 7] # the 2 is omited because we are not using axons
                swc_tree.read_SWC_tree_from_file(filePath, types)
                stats = btmorph.BTStats(swc_tree)
                stats.set_level_to_tree()
                # print ("tree level"+ str(swc_tree.max_height))
                tree_height = swc_tree.max_height
                print str(fileName)+"  "+str(tree_height)
                if tree_height < 0:
                    print(" -> processing " + str(fileName))
                    print("OUT OF RANGE")
                else:
                    try:
                        swc_tree = stats.write_subtree_by_slice_range(slice_begin, slice_end, new_directory,
                                                                      fileName,class_number )
                    except ValueError:
                        print("ERROR  " + str(fileName))

else:
    print ("Enter the required parameters")
