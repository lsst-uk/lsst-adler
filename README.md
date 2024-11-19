# lsst-adler

[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)

[![ci](https://github.com/lsst-uk/lsst-adler/actions/workflows/smoke-test.yml/badge.svg)](https://github.com/lsst-uk/lsst-adler/actions/workflows/smoke-test.yml)
[![pytest](https://github.com/lsst-uk/lsst-adler/actions/workflows/testing-and-coverage.yml/badge.svg)](https://github.com/lsst-uk/lsst-adler/actions/workflows/testing-and-coverage.yml)
[![Read the Docs](https://img.shields.io/readthedocs/adler)](https://adler.readthedocs.io/)
[![PyPI](https://img.shields.io/pypi/v/lsst-adler?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/lsst-adler/)
[![codecov](https://codecov.io/gh/lsst-uk/lsst-adler/branch/main/graph/badge.svg)](https://codecov.io/gh/lsst-uk/lsst-adler)

This project was automatically generated using the LINCC-Frameworks 
[python-project-template](https://github.com/lincc-frameworks/python-project-template).

## Read the Docs

Look here for information on adler, the API and some example notebooks: https://adler.readthedocs.io/en/latest/

## Dev Guide - Getting Started

Before installing any dependencies or writing code, it's a great idea to create a
virtual environment. If you have conda installed locally, you can run the following to
create and activate a new environment.

```
>> conda create --name <env_name> python=3.10
>> conda activate <env_name>
```

Once you have created a new environment, you can install this project for local
development. Git clone the repo and from the top directory run the following commands from within the adler folder:

```
>> pip install -e .'[dev]'
>> pre-commit install
```

**WARNING:** If you're installing on the RSP, then use the following pip command instead:

```
>> pip install --user -e .'[dev]'
```
Furthermore, to get the pre-commit hooks working on RSP, you will probably have to update git within your environment:

```
conda install git
```

If you're also working on the docs:

```
>> conda install pandoc
```

You can then test that everything works by running:

```
adler -s 8268570668335894776
```
This will currently print some phase curve information to the terminal.

One can also read from a local database, for example:

```
adler -s 8268570668335894776 -i tests/data/testing_database.db
```

Notes:
1) The single quotes around `'[dev]'` may not be required for your operating system.
2) `pre-commit install` will initialize pre-commit for this local repository, so
   that a set of tests will be run prior to completing a local commit. For more
   information, see the Python Project Template documentation on 
   [pre-commit](https://lincc-ppt.readthedocs.io/en/latest/practices/precommit.html)
3) Install `pandoc` allows you to verify that automatic rendering of Jupyter notebooks
   into documentation for ReadTheDocs works as expected. For more information, see
   the Python Project Template documentation on
   [Sphinx and Python Notebooks](https://lincc-ppt.readthedocs.io/en/latest/practices/sphinx.html#python-notebooks)

## Dev Guide - Adding notebooks to Read The Docs

- Copy notebook into `docs/notebooks` (N.B. the notebook must have at least one section header* and be using the "Python 3 (ipykernel)" kernel, not some conda env kernel that may only be installed locally)
- Update the toctree in the file `docs/notebooks.rst`
- Ensure necessary requirements are declared in `pyproject.toml` and `docs/requirements.txt`. Also, make sure that the notebook being added to the docs is using the python3 (ipykernel) kernel, not some conda env kernel that may only be installed locally
- To update the docs locally, from the `docs` dir run: `python -m sphinx -T -E -b html -d _build/doctrees -D language=en . ../_readthedocs/html`

* Multiple section headers from a notebook will also show up in the table of contents.

## Dev Guide - Updating pyproject.toml

If you are adding code that requires a new dependency, this needs to be included in pyproject.toml under the `[project]' section:

```
dependencies = [
    "ipykernel", # Support for Jupyter notebooks
    "numpy",
    "lsst-rsp"
    "your-dependency-here"
]
```

If you are adding code that should be run from the command line, this should be set up under `[project.scripts]`:

```
[project.scripts]
adler = "adler.adler:main"
my_command = "adler.module_folder.module_name:function_name"
```

## Dev Guide - Docstrings

The docstrings use the numpydoc format. This is the format expected by the LINCC Frameworks template and results
in docstrings which compile neatly and automatically for the docs. 

An example is below. Don't include the Returns section if your code returns nothing.

```
"""Here is a function that does some cool sciency stuff.

- Perhaps you want to add some bullet points.
- You can do that like this!

Parameters
-----------
arg1 : int
   The first argument to the function.
arg2 : str
   The second argument to the function. Default = "foo".

Returns
----------

return_value : float
   The thing your function returns.

"""

```

For classes:

```
"""A class that contains some important science information.

Attibutes
-----------
attr1 : int
   The first class attribute.

attr2: np.array
   The second class attribute.

"""
```
