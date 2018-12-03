Usage examples
==============

Configuration file
------------------

First you need to define your datasets

.. code-block:: bash

    [cadd]
    type = tabix
    file = cadd_score.tsv.gz
    header = CHR, POS, REF, ALT, RAW, PHRED
    ctypes = str, int, str, str, float, float
    sequence= CHR
    begin = POS
    end = POS

    [exons]
    type = tabix
    file = cds_exons.tsv.gz
    header = CHR, START, STOP, GENE
    ctypes = str, int, int, str
    sequence = CHR
    begin = START
    end = STOP
    indices = GENE,


Group by
--------

You can only use indexes fields as a groupby key. The indexed fields are the ones defined in the **indices** section
of a dataset definition.

::

    gd = Gendas('data/gendas.conf')

    # Total genes
    genes_size = len(set(gd['exons']['GENE']))

    # Group by gene column and compute the mean, maximum and minimum score
    bkg_cadd_mean_by_gene = gd.groupby(gd['exons']['GENE']).aggregate({
        'MEAN': lambda gd: mean(gd['cadd']['PHRED']),
        'MAX': lambda gd: max(gd['cadd']['PHRED']),
        'MIN': lambda gd: min(gd['cadd']['PHRED'])
    })

Merge
-----

You can merge datasets

::

    pan_with_cadd = gd['variants'].merge(gd['cadd'], on=['REF', 'ALT'])







