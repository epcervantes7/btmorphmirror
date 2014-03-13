.. raw:: html

    <style> .red {color:white;background-color:red} </style>
    <style> .green {color:white;background-color:green} </style>


#####################
Validation & testing
#####################

.. _comparison:

Comparison with L-Measure
--------------------------

We compare our library with the "golden standard", L-Measure and the results generated
by L-Measure that are accessible through `NeuroMorpo.org <NeuroMorpho.org>`_.

In most cases the results obtained with the btmorph library are similar; there are some slight differences that reflect slight implementation details and some measures are interpreted differently; implementation details of L-Measure can be found `(here) <http://cng.gmu.edu:8080/Lm/help/index.htm>`_.
We explain the similarities and differences by means of an exemplar analysis performed on one
morphology: `v_e_moto1` `(from here) <http://neuromorpho.org/neuroMorpho/neuron_info.jsp?neuron_name=v_e_moto1>`_. 


.. role:: red
.. role:: green


.. tabularcolumns:: |l|l|p{5cm}|

+---------------------+-----------------+---------------------------+
|Morphometric feature | NeuroMorpho.org | btmorph                   |
+=====================+=================+===========================+
| Soma Surface        | 45216 μm2       | :red:`45238` μm2 [#f0]_   |
+---------------------+-----------------+---------------------------+
| # Stems             | 10              | :green:`10`               |
+---------------------+-----------------+---------------------------+
| # Bifurcations      | 122             | :green:`122`              |
+---------------------+-----------------+---------------------------+
| # Branches          | 254             | :green:`254` [#f1]_       |
+---------------------+-----------------+---------------------------+
| Overall Width       |  1804.67 μm     | 2588.0 μm [#f2]_          |
+---------------------+-----------------+---------------------------+
| Overall Height      |  2259.98 μm     | 2089.0 μm [#f2]_          |
+---------------------+-----------------+---------------------------+
| Overall Depth       |  1701.72 μm     | 2306.0 μm [#f2]_          |
+---------------------+-----------------+---------------------------+
| Average Diameter    |  2.2 μm         | :green:`2.2` μm [#f3]_    |
+---------------------+-----------------+---------------------------+
| Total Length        |  78849.1 μm     | :green:`78849.1` μm       |
+---------------------+-----------------+---------------------------+
| Total Surface       |  512417 μm2     | :green:`512417` μm2       |
+---------------------+-----------------+---------------------------+
| Total Volume        |  390413 μm3     | :green:`390412` μm3       |
+---------------------+-----------------+---------------------------+
| Max Euclidean       |                 |                           |
| Distance            | 765.73 μm       | :red:`1531 μm` [#f4]_     |
+---------------------+-----------------+---------------------------+
| Max Path Distance   | 873.56 μm       | :red:`1817` μm [#f5]_     |
+---------------------+-----------------+---------------------------+
| Max Branch Order    | 3.69            | :green:`3.83` [#f6]_      |
+---------------------+-----------------+---------------------------+
| Average Contraction | 0.94            | :green:`0.9359` [#f7]_    |
+---------------------+-----------------+---------------------------+
| Total Fragmentation | 559             | na [#f8]_                 |
+---------------------+-----------------+---------------------------+
| Partition Asymmetry | 0.43            | :green:`0.43` [#f9]_      |
+---------------------+-----------------+---------------------------+
| Average Rall's      |                 |                           |
| Ratio               |1.25             | :red:`1.69` [#f10]_       |
+---------------------+-----------------+---------------------------+
| Average Bifurcation |                 |                           |
| Angle Local         | 46.83°          | :green:`46.83°`           |
+---------------------+-----------------+---------------------------+
| Average Bifurcation |                 |                           |
| Angle Remote        |  45.74°         | :green:`45.7 °`           |
+---------------------+-----------------+---------------------------+

.. [#f0] In accordance with the three-point soma format, the somatic surface is computed as :math:`A = 4 \times \pi \times r^2`.
.. [#f1] Computed by `stats.no_bifurcations() + stats.no_terminals()`
.. [#f2] We compute the raw, untranslated extend in X,Y and Z dimension. This is different from the translated and truncated extend used by L-Measure.
.. [#f3] Computed by `np.mean(stats.get_diameters())`
.. [#f4] Unclear how the NeuroMorpho.org value is generated. We compute the euclidean distance between each terminal point and the soma. A visual inspection shows that our value is correct.

.. [#f5] See [#f4]_
.. [#f6] This is actually not the maximum as listed on the NeuroMorpho website but the average of either all points, or the bifurcation points.
.. [#f7] Computed as follows: 
:: 

   eds = []
   pls = []
   for node in stats._end_points:
       eds.append(stats.get_segment_Euclidean_length(node))
       pls.append(stats.get_segment_pathlength(node))
   mean(array(eds)/array(pls))

.. [#f8] Not implemented. Depends mostly on the person performing the reconstruction and hence not an important morphometric. Although, this measure can be easily retrieved by counting the compartment in the `get_segment_pathlength` method.

.. [#f9] Computed as follows:
::

   pas = []
   for node in stats._bif_points:
       pas.append(stats.partition_asymmetry(node))
   mean(pas)

.. [#f10] Rall's ratio, :math:`n`, :math:`{D_p}^n={D_{d1}}^n+{D_{d2}}^n` We use two distinct implementations. One based on the L-Measure documentation, namely, brute-force looking for the best :math:`n` over :math:`[0,5]` in 1000 steps. This yielded an average of :math:`n=1.77`, a value that differs from the value reported on NeuroMorpho.org (:math:`n=1.25`). Although, in accordance to the L-Measure documentation many segments are being "discarded" and it is unclear why such a large number is discarded. We only discard the branching points where :math:`d_1 > D_p` or :math:`d_2 > D_p`. In another implementation we use the scipy.optimize module to perform a simplex search for the optimal value of :math:`p` and this yields :math:`p=1.69`. In this second implementation we discard more points, namely points where :math:`p` is not on :math:`[0,5]`.

.. _unit_testing:

Unit testing
------------

Unit-testing refers to testing of elementary pieces of code in a computer program `(Wikipedia) <http://en.wikipedia.org/wiki/Unit_testing>`_. Testing is done using the Python testing framework, called nose tests. In these tests, we compare the outcome of our library to similar outcomes generated by L-Measure that are accessible through the `NeuroMorpho.org <www.neuromorpho.org>`_ website. Note that there are some differences in design and definition of the features as listed :ref:`comparison`.

Unit-tests of this library are provided in the ``tests`` directory and can be run by
::

    nosetests -v tests/stats_test.py

.. note:: Run the unit-tests after change to the code to ensure a) backward compatibility and b) correctness of the results.

