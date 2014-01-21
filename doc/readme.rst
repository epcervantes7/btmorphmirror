btmorph
=======

Introduction
------------

Package to perform analysis of neuronal moprhologies and some associated tricks.

Installation
------------

To run the code, just add the directory containing the files in your PYTHONPATH::

    # in the directory with the *.py files
    export PYTHONPATH=$(pwd):$PYTHONPATH

Test the installation by running the tests:::

    nosetests -v --nocapture tests/tools_test.py
    nosetests -v --nocapture tests/structs_test.py
    nosetests -v --nocapture tests/stats_test.py



Data representation
~~~~~~~~~~~~~~~~~~~

The idea is to represent a morphology ads a tree structure. Tree structures provide an intuitive representation of a morphology and can be easily probed to calculate morphometric features.

The tree is essentially a linked list data structure (STree). Each item in the list/tree is a node (Snode) and contaisn pointers to its parent (``get_parent``) and its children (``get_children``). Each node can store *something*, in this case it should store a location in 3D (``P3D``) that is accessible through the ``SNode.get_content``. Obviously, this tree structure resembles strongly the structure of an SWC file.

Schematically, it looks like this:

.. image:: figures/tree_structure.png
  :scale: 50


Design requirements
-------------------

Quick and basic assessment of morphometric features of neurons. Internal representation based on a tree structure (STree; how to link?). Atomic functions are provided so they can be used in further scripting, analysis and validation. Code is purely serial as single neuron morphometrics go fast anyway. If, however, a batch of morphologies needs to be analyzed, a parallel wrapper can be written.

Input is an SWC file (a filter is provided in btmorphtools to load BBP H5v2 format). The SWC format as used on the curated database NeuroMorpho.org and it is expected that the latest soma-standard is followed, i.e., the soma is described by three points. (see `here <http://neuromorpho.org/neuroMorpho/SomaFormat.html>`_). Tools are provided to filter distinct dendritic types (part of the conversion from BBP H5 format) and hence it can be used to analyze distinct parts of the tree separately when desired.

Morphometrics can be either scalar (= one value per morphology) or vector / distributed (= a distribution of values per morphology). For vector morphometrics, the features can be measures either a branching point, terminal points or both. Other 'points' specified in the SWC file are only used for the internal represention of the geometry.


Morphometric features
~~~~~~~~~~~~~~~~~~~~~

* Scalar: (one per morphological structure under scrutiny)

  * total size: total length of the neurite
  * # stems
  * # branch points
  * # terminal points
  * width (without translation; absolute coordinates; potential extension along the first 3 principal components)
  * height 
  * depth
  * max degree (of neurites sprouting at the soma)
  * max order (of neuritues sprouting at the soma)
  * partition asymmetry (can/cannot be measured at the soma?)

* Vector: (per 'point of interest' PoI):

  * segment path length (incoming)
  * segment euclidean length (incoming)
  * contraction (euclidean / path; incoming)
  * order
  * degree
  * partition asymmetry
  * fractal dimension (of path between soma and PoI)
  * `Clouds`: save x,y,z coordinates for post-hoc histograms analysis or other scalar (e.g., moments) or vector properties (e.g., PCA)


Visualization
~~~~~~~~~~~~~

(simple, using matplotlib):

* Dendrogram
* 2D/3D plot 



Quick example
-------------

::

   import btstructs, btstats, btviz
   import numpy
   import matplotlib.pyplot as plt

   swc_tree = btstructs.STree()
   swc_tree.read_SWC_tree_from_file(file_name)
   stats = btstats.BTStats(swc_tree)

   # get the total length
   total_length = stats.total_length()

   # get the max degree, i.e., degree of the soma
   max_degree = stats.degree_of_node(swc_tree.get_root())

   # generate and save the dendrogram
   btviz.plot_dendrogram(file_name=test_file_name)
   plt.savefig('examplar_dendrogram.pdf')

