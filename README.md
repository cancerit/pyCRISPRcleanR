# pyCRISPRcleanR
| Master                                              | Develop                                               |
| --------------------------------------------------- | ----------------------------------------------------- |
| [![Master Badge][travis-master-badge]][travis-repo] | [![Develop Badge][travis-develop-badge]][travis-repo] |

This is python implementation of Francesco's [CRISPRcleanR] R package for unsupervised identification and
correction of gene independent cell responses to CRISPR-cas9 targeting 

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [Design](#design)
- [Tools](#tools)
	- [pyCRISPRcleanR](#pycrisprcleanr)
	- [inputFormat](#inputformat)
	- [outputFormat](#outputformat)
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

`pyCRISPRcleanR` has multiple commands, listed with `pyCRISPRcleanR --help`.

### pyCRISPRcleanR

Takes the input count data, library file and other associated files/parameters
The output is tab separated files for normalised fold changes and
inverse transformed corrected treatment counts

Various exceptions can occur for malformed input files.

### inputFormat
 * ```gRNA Counts``` file: tab separated file containing following fields
 * sgRNA gene <control_count 1...N> <sample_count 1..N>
 * ```sgRNA library``` file format
 * sgRNA gene chr start end
### outputFormat

  following tab separated output files were produced

 1. crispr_cleanr_normalised_counts.tsv
 * sgRNA: guideRNA
 * gene: gene name as defined in the library file
 * <control sample count:normalised 1..N> : Normalised count
 * <treatment sample count: normalised 1..N> : Normalised count

 2. crispr_cleanr_fold_changes.tsv
 * sgRNA: guideRNA
 * gene: gene name as defined in the library file
 * <treatment sample fold chages: fold changes 1..N>
 * avgFC: average fold change values

 3. crispr_cleanr_corrected_counts.tsv [ generated only when ```--segmentation``` option is selected ]
 * sgRNA: guideRNA
 * gene: gene name as defined in the library file
 * <control sample count:corrected 1..N> : corrected count
 * <treatment sample count:corrected 1..N >: corrected count

 4. crispr_cleanr_alldata.tsv [ generated only when ```--segmentation``` option is selected ]
 * sgRNA: guideRNA
 * <control sample count: raw 1..N> : raw count
 * <treatment sample count: raw 1..N> : raw count
 * gene: gene name as defined in the library file
 * chr: Chromosome name
 * start: gRNA start position
 * end: gRNA end position
 * <control sample count:normalised 1..N> : Normalised count (postfixed _nc)
 * <treatment sample count: normalised 1..N> : Normalised count (postfixed _nc)
 * avgFC: average fold change values
 * BP: Base pair location ( used for DNAcopy analysis)
 * correction: correction factor
 * correctedFC: corrected foldchange values
 * <control sample count:corrected 1..N> : corrected count (postfixed _cc)
 * <treatment sample count:corrected 1..N >: corrected count (postfixed _cc)

## INSTALL
Installing via `pip install`. Simply execute with the path to the compiled 'whl' found on the [release page][pyCRISPRcleanR-releases]:

```bash
pip install pyCRISPRcleanR.X.X.X-py3-none-any.whl
```

Release `.whl` files are generated as part of the release process and can be found on the [release page][pyCRISPRcleanR-releases]

### Package Dependancies

`pip` will install the relevant dependancies, listed here for convenience:
* [NumPy]
* [Pandas]
* [rpy2]
* [plotly]

### R packages

* [DNAcopy] R packages is required to run `pyCRISPRcleanR`.  To facilitate the install process there is
a script `Rsupport/libInstall.R` that can be run to build this for you.

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
pip install --find-links=~/wheels pyCRISPRcleanR
```
<!--refs-->
 [NumPy]: http://www.numpy.org/
 [plotly]: https://plot.ly/python/
 [Pandas]: http://pandas.pydata.org/
 [rpy2]: https://rpy2.bitbucket.io/
 [DNAcopy]: https://www.bioconductor.org/packages/release/bioc/html/DNAcopy.html 
 [CRISPRcleanR]: https://github.com/francescojm/CRISPRcleanR
 [travis-master-badge]: https://travis-ci.org/cancerit/pyCRISPRcleanR.svg?branch=master
 [travis-develop-badge]: https://travis-ci.org/cancerit/pyCRISPRcleanR.svg?branch=develop
 [travis-repo]: https://travis-ci.org/cancerit/pyCRISPRcleanR
 [pyCRISPRcleanR-releases]: https://github.com/cancerit/pyCRISPRcleanR/releases
