"""
File contains:
    P3D2
    SNode2   
    STree2

B. Torben-Nielsen, 2013-10 @ OIST,
2013-01 @ BBP
from BTN legacy code
"""


class P3D2 :
    """
    Basic container to represent and store 3D information
    """
    def __init__(self,xyz,radius,type=7) :
        # """
        # Constructor.

        # Parameters
        # -----------

        # xyz : numpy.array
        #     3D location
        # radius : float
        # type : int
        #     Type asscoiated with the segment according to SWC standards
        # """
        self.xyz = xyz
        self.radius = radius
        self.type = type

    def __str__(self) :
        return "P3DD " + self.xyz + ", R=",self.radius, ", T=",self.type

class SNode2 :
    """
    Simple Node for use with a simple Tree (STree)
    
    By design, the "_content" should be a dictionary. (2013-03-08)
    """
    
    def __init__(self,index) :
        # """
        # Constructor.

        # Parameters
        # -----------
        # index : int
        #    Index, unique name of the SNode2
        # """
        self._parent_node = None
        self._index = index
        self._child_nodes = []
        
    def get_index(self) :
        """
        Return the index of this node

        Returns
        -------
        index : int
        """
        return self._index
        
    def get_parent_node(self) :
        """
        Return the parent node of this one.

        Returns
        -------
        parent : SNode2
           In case of the root, None is returned.Otherwise a SNode2 is returned
        """          
        return self._parent_node
        
    def get_child_nodes(self) :
        """
        Return the child nodes of this one (if any)

        Returns
        -------
        children : list SNode2
           In case of a leaf an empty list is returned
        """                  
        return self._child_nodes
        
    def get_content(self) :
        """
        Return the content dict of a SNode2

        Returns
        -------
        parent : SNode2
           In case of the root, None is returned.Otherwise a SNode2 is returned
        """                  
        return self._content
    
    def set_index(index) :
        """
        Set the unqiue name of a node

        Parameters
        ----------

        index : int
        """
        self._index = index
    
    def set_parent_node(self,parent_node) :
        """
        Set the parent node of a given other node

        Parameters
        ----------
        node : SNode2
        """
        self._parent_node = parent_node
        
    def set_content(self,content) :
        """
        Set the content of a node. The content must be a dict

        Parameters
        ----------
        content : dict
            dict with content. For use in btmorph at least a 'p3d' entry should be present
        """        
        if isinstance(content,dict) :
            self._content = content 
        else :
            raise Exception("SNode2.set_content must receive a dict")    
        
    def add_child(self,child_node) :
        """
        add a child to the children list of a given node

        Parameters
        -----------
        node :  SNode2
        """
        self._child_nodes.append(child_node)
            
    def make_empty(self):
        """
        Clear the node. Unclear why I ever implemented this. Probably to cover up some failed garbage collection
        """
        self._parent_node = None
        self._content = None
        self._child_nodes = []
            
    def remove_child(self, child_node) :
        """
        Remove a child node from the list of children of a specific node

        Parameters
        -----------
        node :  SNode2
            If the child doesn't exist, you get into problems.
        """
        self._child_nodes.remove(child_node)
        
    def __str__(self) :
        return 'SNode2 (ID: '+str(self._index)+')'

    def __lt__(self,other):
        if self._index < other._index :
            return True
    def __le__(self,other):
        if self._index <= other._index :
            return True
    def __gt__(self,other):
        if self._index > other._index :
            return True
    def __ge__(self,other):
        if self._index >= other._index :
            return True
    
    def __copy__(self) : # customization of copy.copy
        ret = SNode2(self._index)
        for child in self._child_nodes :
            ret.add_child(child)
        ret.content = self._content
        ret.set_parent_node(self._parent_node)
        return ret
        
