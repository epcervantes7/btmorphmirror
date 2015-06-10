import glob
import numpy as np

from btmorph import STree2
from btmorph import BTStats

class BTStats_dynamic(object):
    """
    Compute morphometrics over time.
    Input is a "temporal stack" of SWC files
    """
    def __init__(self,partial_name=""):
        """
        Constructor.

        Parameters
        -----------
        partial_name : string
            First part of the file name of the SWC files in the stacks.
            When this string does not include any information abou the
            path, the stack will be sought in the working directory where
            the command is issued.

            .. warning: Make sure the names indicate the order of the files within the stack
        """
        swc_files = sorted(glob.glob(partial_name+"*.swc"))
        self.tree_frames = {}
        count = 0
        for f in swc_files:
            #print("f: {0}".format(f))
            tree = STree2()
            tree.read_SWC_tree_from_file(f)
            self.tree_frames[count]=tree
            count = count +1

    def stem_length(self):
        """
        Compute the stem length over time.

        Returns
        -------
        length_frames : dict
            Dictionary containing a list for each stem. Each list contains
            the length of the stem per time frame.
        """
        ret ={}
        for frame,tree in self.tree_frames.items():
            print(tree)
            children = tree.root.children
            no_stems = 0
            for stem in children:
                if not stem.index in [2,3]: # 2,3 are "fake" compartments part of the 3-point soma specs
                    # recursively find the length 
                    no_stems = no_stems + 1
                    #print("stem: {0}".format(stem))
                    new_children = stem.get_children()
                    
                    #print("new_children: {0}".format(new_children))
                    stem_length = 0
                    p = tree.root
                    while 1 <= len(new_children) < 2:
                        n = new_children[0]
                        stem_length = stem_length + np.sqrt(np.sum((n.content['p3d'].xyz-p.content['p3d'].xyz)**2))
                        new_children = n.get_children()
                        p = n
                        #print("recurse: new_children={0}".format(new_children))
                    #print("Frame {0} has stem_length={1}".format(frame,stem_length))
                    
            #print("Frame {0} has {1} stems".format(frame,no_stems))
        return

    def total_length(self):
        """
        Compute the tota

        Returns
        -------
        length_frames : dict
            Dictionary containing a list for each stem. Each list contains
            the length of the stem per time frame.
        """
        ret ={}
        for frame,tree in self.tree_frames.items():
            stats = BTStats(tree)
            ret[frame]=stats.total_length()
        return ret
                    
        
        
        
