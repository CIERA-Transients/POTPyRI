Installation
============

POTPyRI requires **Python 3.11 or later**.

From PyPI
---------

.. code-block:: bash

    pip install potpyri

From source (development)
--------------------------

1. Clone the repository:

.. code-block:: bash

    git clone https://github.com/CIERA-Transients/POTPyRI
    cd POTPyRI

2. Create and activate a virtual environment (conda or venv). Python 3.12 is recommended.

.. code-block:: bash

    conda create -n potpyri python=3.12
    conda activate potpyri

   # Or with venv:
   # python3.12 -m venv .venv && source .venv/bin/activate  # Linux/macOS

3. Install in editable mode with optional dependencies:

.. code-block:: bash

    pip install -e .              # core package
    pip install -e ".[test]"      # add pytest, pytest-cov for testing
    pip install -e ".[docs]"      # add Sphinx and theme for building docs

The ``-e`` (editable) flag uses the current directory as the live source, so changes are picked up without reinstalling.

Non-Python dependencies
------------------------

The pipeline uses **astrometry.net** and **Source Extractor** (SExtractor). Install them via your system package manager or conda (see the main README). Index files for astrometry.net can be installed with the ``download_anet_index`` script after astrometry.net is on your path.

Building the documentation
---------------------------

With ``.[docs]`` installed, from the project root:

.. code-block:: bash

    cd docs && make html

Output is in ``docs/build/html/``. On Windows use ``make.bat html``.
