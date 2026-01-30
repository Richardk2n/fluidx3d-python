# FluidX3D Python
-----------------
This provides python libraries to work with FluidX3D.
This project was written as a part of my doctoral thesis (will be linked here if it is ever published).
It is used to evaluate the simulation output.

It requires Python 3.8 or newer.

## Installing
-------------

This is not on PyPi or similar.
To install it, clone the repo, enter the directory and run the following command.

```bash
pip install -e .
```

This installs the local project in a way, that allows edits to the local project.

## Developing
-------------

Install development dependencies with (you might want to create a virtualenv first):

```bash
pip install -r requirements.txt
```

The project is formatted with [black](https://github.com/psf/black). You can either configure your IDE to automatically format code with it, run it manually (``black .``) or rely on pre-commit (see below) to format files on git commit.

The project is formatted with [isort](https://github.com/PyCQA/isort). You can either configure your IDE to automatically sort imports with it, run it manually (``isort .``) or rely on pre-commit (see below) to sort files on git commit.

This project uses [pre-commit](https://pre-commit.com/) to enforce code-quality. After cloning the repository install the pre-commit hooks with:

```bash
pre-commit install
```

After that pre-commit will run [all defined hooks](.pre-commit-config.yaml) on every ``git commit`` and keep you from committing if there are any errors.

Contirbutions with incorrect formatting are not accepted.
