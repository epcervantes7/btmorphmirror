"""
File contains:
    P3D2
    SNode2   

B. Torben-Nielsen (legacy code)
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
        #return "P3DD " + str(self.xyz) + ", R=",str(self.radius)+ ", T="+str(self.type)
        return "P3D2 [%.2f %.2f %.2f], R=%.2f" % (self.xyz[0],self.xyz[1],self.xyz[2],self.radius)
        #Return "P3D2 ".format, self.xyz)

class SNode2 :
    """
    Simple Node for use with a simple Tree (STree)
    
    By design, the "content" should be a dictionary. (2013-03-08)
    """
    
    def __init__(self,index) :
        # """
        # Constructor.

        # Parameters
        # -----------
        # index : int
        #    Index, unique name of the SNode2
        # """
        self.parent = None
        self.index = index
        self.children = []
        self.content = {}

    def get_parent(self) :
        """
        Return the parent node of this one.

        Returns
        -------
        parent : :class:`SNode2`
           In case of the root, None is returned.Otherwise a :class:`SNode2` is returned
        """          
        return self.__parent

    def set_parent(self,parent) :
        """
        Set the parent node of a given other node

        Parameters
        ----------
        node : :class:`SNode2`
        """
        self.__parent = parent

    parent = property(get_parent, set_parent)
    
    def get_index(self) :
        """
        Return the index of this node

        Returns
        -------
        index : int
        """
        return self.__index
        
    def set_index(index) :
        """
        Set the unqiue name of a node

        Parameters
        ----------

        index : int
        """
        self.__index = index
    
    index = property(get_index, set_index)

    def get_children(self) :
        """
        Return the children nodes of this one (if any)

        Returns
        -------
        children : list :class:`SNode2`
           In case of a leaf an empty list is returned
        """                  
        return self.__children
        
    def set_children(self, children) :
        """
        Set the children nodes of this one

        Parameters
        ----------

        children: list :class:`SNode2`
        """
        self.__children = children

    children = property(get_children, set_children)

    def get_content(self) :
        """
        Return the content dict of a :class:`SNode2`

        Returns
        -------
        parent : :class:`SNode2`
           In case of the root, None is returned.Otherwise a :class:`SNode2` is returned
        """                  
        return self.__content
    
    def set_content(self,content) :
        """
        Set the content of a node. The content must be a dict

        Parameters
        ----------
        content : dict
            dict with content. For use in btmorph at least a 'p3d' entry should be present
        """        
        if isinstance(content,dict) :
            self.__content = content 
        else :
            raise Exception("SNode2.set_content must receive a dict")    

    content = property(get_content, set_content)

    def add_child(self,child_node) :
        """
        add a child to the children list of a given node

        Parameters
        -----------
        node :  :class:`SNode2`
        """
        self.children.append(child_node)
            
    def make_empty(self):
        """
        Clear the node. Unclear why I ever implemented this. Probably to cover up some failed garbage collection
        """
        self.parent = None
        self.content = {}
        self.children = []
            
    def remove_child(self, child) :
        """
        Remove a child node from the list of children of a specific node

        Parameters
        -----------
        node :  :class:`SNode2`
            If the child doesn't exist, you get into problems.
        """
        self.children.remove(child)
        
    def __str__(self) :
        return 'SNode2 (ID: '+str(self.index)+')'

    def __lt__(self,other):
        if self.index < other.index :
            return True
    def __le__(self,other):
        if self.index <= other.index :
            return True
    def __gt__(self,other):
        if self.index > other.index :
            return True
    def __ge__(self,other):
        if self.index >= other.index :
            return True
    
    def __copy__(self) : # customization of copy.copy
        ret = SNode2(self.index)
        for child in self.children :
            ret.add_child(child)
        ret.content = self.content
        ret.parent = self.parent
        return ret

