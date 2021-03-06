.. note::
    :class: sphx-glr-download-link-note

    Click :ref:`here <sphx_glr_download_auto_examples_plotting_plot_celltype_quantification.py>` to download the full example code
.. rst-class:: sphx-glr-example-title

.. _sphx_glr_auto_examples_plotting_plot_celltype_quantification.py:


visualize cell fractions
========================

This example demonstrates how to generate celltype quantification plots. These types of plots 
can be used to visually represent the number of cells that belong to a certain subset or condition.




.. code-block:: python


    import besca as bc 

    #import dataset to workwith
    adata = bc.datasets.pbmc_storage_processed_bbknn()







quantify specific celllabels as a stacked barplot



.. code-block:: python


    bc.pl.celllabel_quant_stackedbar(adata, count_variable = 'donor', subset_variable = 'storage_condition');




.. image:: /auto_examples/plotting/images/sphx_glr_plot_celltype_quantification_001.png
    :class: sphx-glr-single-img




quantify number of cells in each louvain cluster belonging to a specific subset variable (e.g. donor)

Note that the louvain clusters are brought into the correct order.



.. code-block:: python


    bc.pl.louvain_quant_stackedbar(adata, subset_variable = 'donor');




.. image:: /auto_examples/plotting/images/sphx_glr_plot_celltype_quantification_002.png
    :class: sphx-glr-single-img




quantify number of cells belong to each condition in a specific subset

here each dot represents one donor, the boxplots are grouped according to storage condition



.. code-block:: python


    bc.pl.celllabel_quant_boxplot(adata, count_variable = 'louvain', subset_variable = 'donor', condition_identifier = 'storage_condition',  plot_percentage = True);




.. image:: /auto_examples/plotting/images/sphx_glr_plot_celltype_quantification_003.png
    :class: sphx-glr-single-img




here you can also choose to plot total counts instead of percentages



.. code-block:: python


    bc.pl.celllabel_quant_boxplot(adata, count_variable = 'louvain', subset_variable = 'donor', condition_identifier = 'storage_condition',  plot_percentage = False);



.. image:: /auto_examples/plotting/images/sphx_glr_plot_celltype_quantification_004.png
    :class: sphx-glr-single-img




**Total running time of the script:** ( 0 minutes  5.903 seconds)


.. _sphx_glr_download_auto_examples_plotting_plot_celltype_quantification.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download

     :download:`Download Python source code: plot_celltype_quantification.py <plot_celltype_quantification.py>`



  .. container:: sphx-glr-download

     :download:`Download Jupyter notebook: plot_celltype_quantification.ipynb <plot_celltype_quantification.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.readthedocs.io>`_
