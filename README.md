# pyCRISPRCleanR
This is python implementation of Francesco's [CRISPRclearR] R package for unsupervised identification and
correction of gene independent cell responses to CRISPR-cas9 targeting 

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Design](#design)
- [Tools](#tools)
	- [pyCRISPRCleanR](#pycrisprcleanr)
- [INSTALL](#install)
	- [Package Dependencies](#package-dependencies)
  - [R packages](#r-packages)
- [Development environment](#development-environment)
	- [Development Dependencies](#development-dependencies)
		- [Setup VirtualEnv](#setup-virtualenv)
- [Cutting a release](#cutting-a-release)
	- [Install via `.whl` (wheel)](#install-via-whl-wheel)

<!-- /TOC -->

## Design
Uses DNAcopy R pcakage to perform CBS[ Circular Binary Segmentation of count  data ]

## Tools

`pyCRISPRCleanR` has multiple commands, listed with `pyCRISPRCleanR --help`.

### pyCRISPRCleanR

Takes the input count data, library file and other associated files/parameters
The output is tab separated files for normalised fold changes and
inverse transformed corrected treatment counts

Various exceptions can occur for malformed input files.

## INSTALL
Installing via `pip install` .Simply execute with the path to the compiled 'whl' found on the [release page][pyCRISPRCleanR-releases]:

```bash
pip install pyCRISPRCleanR.X.X.X-py3-none-any.whl
```

Release `.whl` files are generated as part of the release process and can be found on the [release page][pyCRISPRCleanR-releases]

### Package Dependancies

`pip` will install the relevant dependancies, listed here for convenience:
* [NumPy]
* [Pandas]
* [rpy2]

### R packages

* [DNAcopy] R packages is required to run `pyCRISPRcleanR`.  To facilitate the install process there is
a script `Rsupport/libInstall.R` that can be run to build this for you.  See this file for the list
of packages.

Alternatively you can run:

```
cd Rsupport
./setupR.sh path_to_install_to
```

Appending `1` to the command to request a complete local build of `R` (3.3.0).

## Development environment

This project uses git pre-commit hooks.  As these will execute on your system it
is entirely up to you if you activate them.

If you want tests, coverage reports and lint-ing to automatically execute before
a commit you can activate them by running:

```
git config core.hooksPath git-hooks
```

Only a test failure will block a commit, lint-ing is not enforced (but please consider
following the guidance).

You can run the same checks manually without a commit by executing the following
in the base of the clone:

```bash
./run_tests.sh
```

### Development Dependencies

#### Setup VirtualEnv

```
cd $PROJECTROOT
hash virtualenv || pip3 install virtualenv
virtualenv -p python3 env
source env/bin/activate
python setup.py develop # so bin scripts can find module
```

For testing/coverage (`./run_tests.sh`)

```
source env/bin/activate # if not already in env
pip install pytest
pip install pytest-cov
```

__Also see__ [Package Dependancies](#package-dependancies)

### Cutting a release

__Make sure the version is incremented__ in `./setup.py`

### Install via `.whl` (wheel)

Generate `.whl`

```bash
source env/bin/activate # if not already
python setup.py bdist_wheel -d dist
```

Install .whl

```bash
# this creates an wheel archive which can be copied to a deployment location, e.g.
scp dist/pyCRISPRCleanR.X.X.X-py3-none-any.whl user@host:~/wheels
# on host
pip install --find-links=~/wheels pyCRISPRCleanR
```
<!--refs-->
 [NumPy]: http://www.numpy.org/
 [Pandas]: http://pandas.pydata.org/
 [rpy2]: https://rpy2.bitbucket.io/
 [DNAcopy]: https://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html 
