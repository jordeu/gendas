Gendas
======

Gendas is a "pandas" like library focus on genomic datasets. It's a data query engine highly focus to work with
genomic dataset that they can always be map to a genomic position and a sequence id.

The two main goals that differentiate Gendas from other generic data query engines are:

- **Native parallelization:** When we apply some data manipulations in a genomic dataset, most of the times
  it's possible to split the genome in segments and compute each segment in parallel. Gendas tries to hide
  the parallelization complexity to the end users using a declarative syntax (very similar to Pandas and
  other dataframe libraries).


- **Native data sources:** In the genomics field there are several widely use data formats (tabix, bams, fasta,
  fastq, bed...) Gendas can directly work and mix this data formats without reading everything in memory
  or inserting all the data to a different data storage system.


Installation
------------

Install the latest release using pip::

        pip install https://github.com/jordeu/gendas/archive/0.1.tar.gz

Or the development version::

        pip install git+https://github.com/jordeu/gendas


.. note::

    Gendas is still in a very early alpha development stage. Please use it as a proof of concept and to
    help use debug the current implementation and to improve the interaface.


Usage
-----

Check our `example scripts and notebooks <examples>`_ to learn how to use it. And browse the documentation at `readthedocs <https://gendas.readthedocs.io/en/latest/introduction.html>`_.


License
-------

Gendas is available to the general public subject to certain conditions described in its `LICENSE <LICENSE.txt>`_.
