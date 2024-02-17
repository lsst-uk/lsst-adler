# adler

[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)

[![ci](https://github.com/lsst-uk/adler/actions/workflows/smoke-test.yml/badge.svg)](https://github.com/lsst-uk/adler/actions/workflows/smoke-test.yml)
[![pytest](https://github.com/lsst-uk/adler/actions/workflows/testing-and-coverage.yml/badge.svg)](https://github.com/lsst-uk/adler/actions/workflows/testing-and-coverage.yml)
[![Read the Docs](https://img.shields.io/readthedocs/adler)](https://adler.readthedocs.io/)
[![PyPI](https://img.shields.io/pypi/v/adler?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/adler/)
[![codecov](https://codecov.io/gh/lsst-uk/adler/branch/main/graph/badge.svg)](https://codecov.io/gh/lsst-uk/adler)

This project was automatically generated using the LINCC-Frameworks 
[python-project-template](https://github.com/lincc-frameworks/python-project-template).


## Dev Guide - Getting Started

Before installing any dependencies or writing code, it's a great idea to create a
virtual environment. If you have conda installed locally, you can run the following to
create and activate a new environment.

```
>> conda create env -n <env_name> python=3.10
>> conda activate <env_name>
```

Once you have created a new environment, you can install this project for local
development using the following commands from within the adler folder:

```
>> pip install -e .'[dev]'
>> pre-commit install
```

**WARNING:** If you're installing on the RSP, then use the following pip command instead:

```
>> pip install --user -e .'[dev]'
```

If you're also working on the docs:

```
>> conda install pandoc
```

You can then test that everything works by running:

```
adler -s 8268570668335894776
```
This currently prints a friendly message to the terminal.

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
