Installation
============

1. Clone the repository:

.. code-block:: bash

    git clone https://github.com/CIERA-Transients/POTPyRI

2. Create and activate a virtual environment. Anaconda and Python 3.12 are recommended, but not required.

.. code-block:: bash

    conda create -n "potpyri" python=3.12
    conda activate potpyri

3. Install live version with pip:

.. code-block:: bash

   cd POTPyRI
   pip install -e .

Specifically, the ``-e`` (editable) flag uses the specified directory as the "live" source code, and does not copy it to the usual site-packages directory.
