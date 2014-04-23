"""
Basic visualization of neurite morphologies using matplotlib.

Usage is restricted to morphologies in the sWC format with the three-point soma `standard <http://neuromorpho.org/neuroMorpho/SomaFormat.html>`_

B. Torben-Nielsen
"""
import sys,time
sys.setrecursionlimit(10000)

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec


import btmorph

""" internal constants required for the dendrogram generation """
H_SPACE = 20
V_SPACE = 0
C = 'k'

max_width = 0
max_height = 0

def _plot_3D_figure(SWC,my_color_list = ['r','g','b','c','m','y','r--','b--','g--']):
    """
    Parameters
    -----------
    SWC: dict
    """
    for index in SWC.keys() : # not ordered but that has little importance here
        # draw a line segment from parent to current point
        current_SWC = SWC[index]
        #print 'index: ', index, ' -> ', current_SWC
        c_x = current_SWC[0]
        c_y = current_SWC[1]
        c_z = current_SWC[2]
        c_r = current_SWC[3]
        parent_index = current_SWC[4]
                
        if(index <= 3) :
            #print 'do not draw the soma and its CNG, !!! 2 !!! point descriptions'
            pass
        else :
            parent_SWC = SWC[parent_index]
            p_x = parent_SWC[0]
            p_y = parent_SWC[1]
            p_z = parent_SWC[2]
            p_r = parent_SWC[3]
            # print 'index:', index, ', len(cs)=', len(cs)
            pl = plt.plot([p_x,c_x],[p_y,c_y],[p_z,c_z],my_color_list[current_SWC[5]-1],linewidth=c_r)

def true_2D_projections(file_name='P20-DEV139.CNG.swc',align=True,outN=None,bar=None,depth="Z"):
    """
    Three 2D projections

    Parameters
    -----------
    file_name : string
        File name of the SWC file to plots
    depth : string
        Set which axis represents "depth". In experimental work, the \
        Z-axis is depth (as in my PPNeurMorphC) but in NeuroMorpho the \
        Y-axis is the depth. (Depth is used to denote the axis from deep to superficial)
    align : boolean
        Translate the figure so that the soma is on the origin [0,0,0].
    outN : string
        File name of the output file. Extension of this file sets the file type
    bar : list of int or real
        Three values to set the thicks and marks on the plot in X,Y and Z-dimension
    depth : string
        Set the axis representing the depth (= axis from superficial to deep). In most SWC files \
        this is 'Y'. The other option is 'Z', that is more consistent with the usual Cartesian coordinate systems

    """
    my_color_list = ['r','g','b','c','m','y','k','g','DarkGray']
    # read the SWC into a dictionary: key=index, value=(x,y,z,d,parent)
    x = open(file_name,'r')
    soma = None
    SWC = {}
    for line in x :
        if(not line.startswith('#')) :
            splits = line.split()
            index = int(splits[0])
            if index == 1:
                soma_x = float(splits[2])
                soma_y = float(splits[3])
                soma_z = float(splits[4])

            n_type = int(splits[1])
            if not align:
                x = float(splits[2])
                y = float(splits[3])
                z = float(splits[4])
            else:
                x = float(splits[2])-soma_x
                y = float(splits[3])-soma_y
                z = float(splits[4])-soma_z         
            r = float(splits[5])
            parent = int(splits[-1])
            SWC[index] = (x,y,z,r,parent,n_type)

    fig = plt.figure()
    
    if depth == "Y":
        ax = plt.subplot2grid((2,2), (0, 0))
        for index in SWC.keys():
            if index < 3:
                continue
            C = SWC[index]
            P = SWC[C[4]] # C[4] is parent index
            pl = plt.plot([P[0],C[0]],[P[1],C[1]],my_color_list[C[5]-1],linewidth=C[3])
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.xticks([0,bar[0]])
        plt.yticks([0,bar[1]])

        ax2 = plt.subplot2grid((2,2), (0, 1))
        for index in SWC.keys():
            if index < 3:
                continue
            C = SWC[index]
            P = SWC[C[4]] # C[4] is parent index
            pl2 = plt.plot([P[2],C[2]],[P[1],C[1]],my_color_list[C[5]-1],linewidth=C[3])
        plt.xlabel("Z")
        plt.ylabel("Y")
        plt.xticks([0,bar[2]])
        plt.yticks([0,bar[1]])
        
        ax3 = plt.subplot2grid((2,2), (1, 1))
        for index in SWC.keys():
            if index < 3:
                continue
            C = SWC[index]
            P = SWC[C[4]] # C[4] is parent index            
            pl3 = plt.plot([P[0],C[0]],[P[2],C[2]],my_color_list[C[5]-1],linewidth=C[3])
        plt.xlabel("X")
        plt.ylabel("Z")
        plt.xticks([0,bar[0]])
        plt.yticks([0,bar[2]])
    else: # default for my code: Z is depth
        ax = plt.subplot2grid((2,2), (0, 0))
        for index in SWC.keys():
            if index < 3:
                continue
            C = SWC[index]
            P = SWC[C[4]] # C[4] is parent index
            pl = plt.plot([P[0],C[0]],[P[2],C[2]],my_color_list[C[5]-1],linewidth=C[3])
        plt.xlabel("X")
        plt.ylabel("Z")
        plt.xticks([0,bar[0]])
        plt.yticks([0,bar[2]])

        ax2 = plt.subplot2grid((2,2), (0, 1))
        for index in SWC.keys():
            if index < 3:
                continue
            C = SWC[index]
            P = SWC[C[4]] # C[4] is parent index
            pl2 = plt.plot([P[1],C[1]],[P[2],C[2]],my_color_list[C[5]-1],linewidth=C[3])
        plt.xlabel("Y")
        plt.ylabel("Z")
        plt.xticks([0,bar[1]])
        plt.yticks([0,bar[2]])
        
        ax3 = plt.subplot2grid((2,2), (1, 1))
        for index in SWC.keys():
            if index < 3:
                continue
            C = SWC[index]
            P = SWC[C[4]] # C[4] is parent index            
            pl3 = plt.plot([P[0],C[0]],[P[1],C[1]],my_color_list[C[5]-1],linewidth=C[3])
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.xticks([0,bar[0]])
        plt.yticks([0,bar[1]])

    plt.tight_layout()

    if not outN == None:
        plt.savefig(outN)

def plot_2D_SWC(file_name='P20-DEV139.CNG.swc',cs=None,synapses=None,syn_cs=None,outN=None,align=True,offset=None,show_axis=False,XZ=False,filter=range(10)) :
    """
    2D matplotlib plot of a neuronal moprhology. Projection can be in XY and XZ.
    The SWC has to be formatted with a "three point soma".
    Colors can be provided

    Parameters
    -----------
    file_name : string
        File name of the SWC file to plots
    cs : list of float
        Raw values that will be plotted on the morphology according to a colormap
    synapses : list of int
        Indices of the compartments where synapses (= small circles) should be drawn
    syn_c : string
        Color of the synapses ('r','g', 'b', ...)
    outN : string
        File name of the output file. Extension of this file sets the file type
    align : boolean
        Translate the figure so that the soma is on the origin [0,0,0].
    offset : list on float
        List of length 3 with X,Y and Z shift of the morphology to be plotted. Not to be used in conjunction with the "align" option
    show_axis : boolean
        whether or not to draw the axis
    filter : list
        List of integers indicating the SWC types to be included (1:soma, 2:axon, 3:basal,4:apical,...). By default, all SWC types are included
    XZ : boolean
        Default False means the XY projection is drawn. If True, the XZ projection is drawn

    """
    # resolve some potentially conflicting arguments
    if not offset == None:
        align = False
    # read the SWC into a dictionary: key=index, value=(x,y,z,d,parent)
    x = open(file_name,'r')
    SWC = {}
    for line in x :
        if(not line.startswith('#')) :
            splits = line.split()
            index = int(splits[0])
            if index == 1:
                soma_x = float(splits[2])
                soma_y = float(splits[3])
                soma_z = float(splits[4])
            
            n_type = int(splits[1])
            if not offset == None :
                x = float(splits[2])-offset[0]
                y = float(splits[3])-offset[1]
                z = float(splits[4])-offset[2]
            elif not align:
                x = float(splits[2])
                y = float(splits[3])
                z = float(splits[4])
            else:
                x = float(splits[2])-soma_x
                y = float(splits[3])-soma_y
                z = float(splits[4])-soma_z 
            r = float(splits[5])
            parent = int(splits[-1])
            if n_type in filter:
                SWC[index] = (x,y,z,r,parent,n_type)

    my_color_list = ['r','g','b','c','m','y','r--','b--','g--']
            
    if cs == None :
        pass
    else :
        norm = colors.normalize(np.min(cs),np.max(cs))
        Z = [[0,0],[0,0]]
        levels=range(int(np.min(cs)),int(np.max(cs)),100)
        levels = np.linspace(np.min(cs),np.max(cs),100)
        CS3 = plt.contourf(Z,levels,cmap=cm.jet)
        plt.clf()

    min_y = 100000.0

    syn_index = -1
    for index in SWC.keys() : # not ordered but that has little importance here
        # draw a line segment from parent to current point
        current_SWC = SWC[index]
        #print 'index: ', index, ' -> ', current_SWC
        c_x = current_SWC[0]
        c_y = current_SWC[1]
        c_z = current_SWC[2]
        c_r = current_SWC[3]
        parent_index = current_SWC[4]

        if(c_y < min_y) :
            min_y = c_y
                
        if(index <= 3) :
            #print 'do not draw the soma and its CNG, 2 point descriptions'
            pass
        else :
            parent_SWC = SWC[parent_index]
            p_x = parent_SWC[0]
            p_y = parent_SWC[1]
            p_z = parent_SWC[2]
            p_r = parent_SWC[3]
            if(p_y < min_y) :
                min_y= p_y
            # print 'index:', index, ', len(cs)=', len(cs)
            if(cs == None) :
                if XZ:
                    pl = plt.plot([p_z,c_z],[p_x,c_x],my_color_list[current_SWC[5]-1],linewidth=c_r/2.0)
                    #pl = plt.plot([p_x,c_x],[p_z,c_z],my_color_list[current_SWC[5]-1],linewidth=c_r)
                else:
                    pl = plt.plot([p_x,c_x],[p_y,c_y],my_color_list[current_SWC[5]-1],linewidth=c_r/2.0)
            else :
                try :
                    if XZ:
                        pl = plt.plot([p_x,c_x],[p_z,c_z],my_color_list[current_SWC[5]-1],linewidth=c_r/2.0)
                    else:
                        pl = plt.plot([p_x,c_x],[p_y,c_y],my_color_list[current_SWC[5]-1],linewidth=c_r/2.0)                    
                except Exception :
                    pass# it's ok: it's the list size...

        # add the synapses
        if(synapses != None) :
            if index in synapses :
                syn_index = syn_index + 1
                if syn_cs :
                    if XZ:
                        plt.plot(c_x,c_z,syn_cs)
                    else:
                        plt.plot(c_x,c_y,syn_cs)
                    
                else :
                    if XZ:
                        plt.plot(c_x,c_z,syn_cs,'ro')
                    else:
                        plt.plot(c_x,c_y,syn_cs,'ro')
                    


    frame1 = plt.gca()

    if show_axis :
        frame1.axes.get_xaxis().set_ticks([])
        frame1.axes.get_yaxis().set_ticks([])
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        
    # # draw a scale bar
    # if offset == None :
    #     plt.plot([0,100],[min_y*1.1,min_y*1.1],'k',linewidth=5) # 250 for MN, 100 for granule
    # else :
    #     pass # when offset is used, the figure is to scale and no scalebar needed
    
    if(cs != None) :
        cb = plt.colorbar(CS3) # bit of a workaround, but it seems to work
        ticks_f = np.linspace(np.min(cs)-1,np.max(cs)+1,5)
        ticks_i = map(int,ticks_f)
        cb.set_ticks(ticks_i)

    if(outN != None) :
        plt.savefig(outN)


def plot_3D_SWC(file_name='P20-DEV139.CNG.swc',cs=None,synapses=None,syn_cs=None,outN=None,offset=None,align=True,filter=range(10)) :
    """
    3D matplotlib plot of a neuronal morphology. The SWC has to be formatted with a "three point soma".
    Colors can be provided and synapse location marked

    Parameters
    -----------
    file_name : string
        File name of the SWC file to plots
    cs : list of float
        Raw values that will be plotted on the morphology according to a colormap
    synapses : list of int
        Indices of the compartments where synapses (= small circles) should be drawn
    syn_cs : string
        Color of the synapses ('r','g', 'b', ...)
    outN : string
        File name of the output file. Extension of this file sets the file type
    offset : list on float
        List of length 3 with X,Y and Z shift of the morphology to be plotted. Not to be used in conjunction with the "align" option
    show_axis : boolean
        whether or not to draw the axis
    filter : list
        List of integers indicating the SWC types to be included (1:soma, 2:axon, 3:basal,4:apical,...). By default, all SWC types are included        

    """
    my_color_list = ['r','g','b','c','m','y','r--','b--','g--']

    # resolve some potentially conflicting arguments
    if not offset == None:
        align = False    
    
    if(cs == None) :
        pass
    else :
        norm = colors.normalize(np.min(cs),np.max(cs))
        Z = [[0,0],[0,0]]
        levels=range(int(np.min(cs)),int(np.max(cs)),100)
        levels = np.linspace(np.min(cs),np.max(cs),100)
        CS3 = plt.contourf(Z,levels,cmap=cm.jet)
        plt.clf()
    
    # read the SWC into a dictionary: key=index, value=(x,y,z,d,parent)
    x = open(file_name,'r')
    SWC = {}
    # for line in x :
    #     if(not line.startswith('#')) :
    #         splits = line.split()
    #         index = int(splits[0])
    #         n_type = int(splits[1])
    #         x = float(splits[2])
    #         y = float(splits[3])
    #         z = float(splits[4])
    #         r = float(splits[5])
    #         parent = int(splits[-1])
    #         SWC[index] = (x,y,z,r,parent,n_type)

    for line in x :
        if(not line.startswith('#')) :
            splits = line.split()
            index = int(splits[0])
            if index == 1:
                soma_x = float(splits[2])
                soma_y = float(splits[3])
                soma_z = float(splits[4])
            
            n_type = int(splits[1])
            if not offset == None :
                x = float(splits[2])-offset[0]
                y = float(splits[3])-offset[1]
                z = float(splits[4])-offset[2]
            elif not align:
                x = float(splits[2])
                y = float(splits[3])
                z = float(splits[4])
            else:
                x = float(splits[2])-soma_x
                y = float(splits[3])-soma_y
                z = float(splits[4])-soma_z 
            r = float(splits[5])
            parent = int(splits[-1])
            if n_type in filter:
                SWC[index] = (x,y,z,r,parent,n_type)                        

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    for index in SWC.keys() : # not ordered but that has little importance here
        # draw a line segment from parent to current point
        current_SWC = SWC[index]
        #print 'index: ', index, ' -> ', current_SWC
        c_x = current_SWC[0]
        c_y = current_SWC[1]
        c_z = current_SWC[2]
        c_r = current_SWC[3]
        parent_index = current_SWC[4]
                
        if(index <= 3) :
            print 'do not draw the soma and its CNG, !!! 2 !!! point descriptions'
        else :
            parent_SWC = SWC[parent_index]
            p_x = parent_SWC[0]
            p_y = parent_SWC[1]
            p_z = parent_SWC[2]
            p_r = parent_SWC[3]
            # print 'index:', index, ', len(cs)=', len(cs)
            if cs == None :
                pl = plt.plot([p_x,c_x],[p_y,c_y],[p_z,c_z],my_color_list[current_SWC[5]-1],linewidth=c_r/2.0)
            else :
                try :
                    pl = plt.plot([p_x,c_x],[p_y,c_y], \
                                  c=cm.jet(norm(cs[index])),linewidth=c_r)
                except Exception :
                    print 'something going wrong here'
                    # pass# it's ok: it's the list size...

        # add the synapses
        if(synapses != None) :
            if(index in synapses) :
                if syn_cs :
                    plt.plot(c_x,c_y,syn_cs)
                else :
                    plt.plot(c_x,c_y,'ro')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

    if(cs != None) :
        cb = plt.colorbar(CS3) # bit of a workaround, but it seems to work
        ticks_f = np.linspace(np.min(cs)-1,np.max(cs)+1,5)
        ticks_i = map(int,ticks_f)
        cb.set_ticks(ticks_i)

    if(outN != None) :
        plt.savefig(outN)

    return ax

def plot_dendrogram(file_name,transform='plain',shift=0,c='k',radius=True,rm=20000.0,ra=200,outN=None) :
    global C, RM, RA, max_width, max_height # n.a.s.t.y.
    '''
    Generate a dendrogram from an SWC file. The SWC has to be formatted with a "three point soma"

    Parameters
    -----------
    file_name : string
        File name of the SWC file to plots
    transform : string
        Either 'plain' or 'lambda'. Plain means no transform while 'lambda' performs an elecotrtonic transform
    shift : float
        Offset in the x-direction
    c : string
        Color ('r','g', 'b', ...)
    radius : boolean
        Plot a wire (False) dendrogram or one with the thickness of the processes (True)
    rm : float
       Membrane resistance. Only needed when transform = 'lambda'
    rm : float
       Axial resistance. Only needed when transform = 'lambda'
    outN : string
        File name of the output file. Extension of this file sets the file type

    '''
    swc_tree = btmorph.STree2()
    swc_tree = swc_tree.read_SWC_tree_from_file(file_name)
    RM = rm
    RA = ra
    C = c
    max_height = 0
    max_width = 0
    plt.clf()
    print 'Going to build the dendrogram. This might take a while...'
    ttt = time.time()
    _expand_dendrogram(swc_tree.get_root(),swc_tree,shift,0,radius=radius,transform=transform)
    if(transform == 'plain') :
        plt.ylabel('L (micron)')
    elif(transform == 'lambda') :
        plt.ylabel('L (lambda)')
    print (time.time() - ttt), ' later the dendrogram was finished. '

    print 'max_widht=%f, max_height=%f' % (max_width,max_height)
    x_bound = (max_width / 2.0) + (0.1*max_width)
    max_y_bound = max_height + 0.1*max_height
    plt.axis([-1.0*x_bound,x_bound,-0.1*max_height,max_y_bound])

    plt.plot([x_bound,x_bound],[0,100],'k', linewidth=5) # 250 for MN, 100 for granule

    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_visible(False)
    frame1.axes.get_yaxis().set_visible(False)
    
    if(outN != None) :
        plt.savefig(outN)

def _expand_dendrogram(cNode,swc_tree,off_x,off_y,radius,transform='plain') :
    global max_width,max_height # middle name d.i.r.t.y.
    '''
    Gold old fashioned recursion... sys.setrecursionlimit()!
    '''
    place_holder_h = H_SPACE
    max_degree = swc_tree.degree_of_node(cNode)
    required_h_space = max_degree * place_holder_h
    start_x = off_x-(required_h_space/2.0)
    if(required_h_space > max_width) :
        max_width = required_h_space

    
    children = cNode.get_child_nodes()

    if swc_tree.is_root(cNode) :
        print 'i am expanding the root'
        children.remove(swc_tree.get_node_with_index(2))
        children.remove(swc_tree.get_node_with_index(3))
    
    for cChild in  children :
        l = _path_between(swc_tree,cChild,cNode,transform=transform)
        r = cChild.get_content()['p3d'].radius

        cChild_degree = swc_tree.degree_of_node(cChild)
        new_off_x = start_x + ( (cChild_degree/2.0)*place_holder_h )
        new_off_y = off_y+(V_SPACE*2)+l
        r = r if radius  else 1
        plt.vlines(new_off_x,off_y+V_SPACE,new_off_y,linewidth=r,colors=C)
        if((off_y+(V_SPACE*2)+l) > max_height) :
            max_height = off_y+(V_SPACE*2)+l

        _expand_dendrogram(cChild,swc_tree,new_off_x,new_off_y,radius=radius,transform=transform)

        start_x = start_x + (cChild_degree*place_holder_h)
        plt.hlines(off_y+V_SPACE,off_x,new_off_x,colors=C)

def _path_between(swc_tree,deep,high,transform='plain') :
    path = swc_tree.path_to_root(deep)
    pl = 0
    pNode = deep
    for node in path[1:] :
        pPos = pNode.get_content()['p3d'].xyz
        cPos = node.get_content()['p3d'].xyz
        pl = pl + np.sqrt(np.sum((cPos-pPos)**2))
        #pl += np.sqrt( (pPos.x - cPos.x)**2 + (pPos.y - cPos.y)**2 + (pPos.z - cPos.z)**2 )
        pNode = node
        if(node == high) : break
        
    if(transform == 'plain'):
        return pl
    elif(transform == 'lambda') :
        DIAM = (deep.get_content()['p3d'].radius*2.0 + high.get_content()['p3d'].radius*2.0) /2.0 # naive...
        c_lambda = np.sqrt(1e+4*(DIAM/4.0)*(RM/RA))
        return pl / c_lambda
