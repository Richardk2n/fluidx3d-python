FluidX3D Python
======================
#.. image:: https://badge.fury.io/py/pylsp-mypy.svg
#    :target: https://badge.fury.io/py/pylsp-mypy

#.. image:: https://github.com/python-lsp/pylsp-mypy/workflows/Python%20package/badge.svg?branch=master
#    :target: https://github.com/python-lsp/pylsp-mypy/

This provides python libraries to work with FluidX3D

It requires Python 3.8 or newer.


Developing
-------------

Install development dependencies with (you might want to create a virtualenv first):

::

   pip install -r requirements.txt

The project is formatted with `black`_. You can either configure your IDE to automatically format code with it, run it manually (``black .``) or rely on pre-commit (see below) to format files on git commit.

The project is formatted with `isort`_. You can either configure your IDE to automatically sort imports with it, run it manually (``isort .``) or rely on pre-commit (see below) to sort files on git commit.

The project uses two rst tests in order to assure uploadability to pypi: `rst-linter`_ as a pre-commit hook and `rstcheck`_ in a GitHub workflow. This does not catch all errors.

This project uses `pre-commit`_ to enforce code-quality. After cloning the repository install the pre-commit hooks with:

::

   pre-commit install

After that pre-commit will run `all defined hooks`_ on every ``git commit`` and keep you from committing if there are any errors.

.. _black: https://github.com/psf/black
.. _isort: https://github.com/PyCQA/isort
.. _rst-linter: https://github.com/Lucas-C/pre-commit-hooks-markup
.. _rstcheck: https://github.com/myint/rstcheck
.. _pre-commit: https://pre-commit.com/
.. _all defined hooks: .pre-commit-config.yaml
