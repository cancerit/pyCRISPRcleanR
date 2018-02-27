# archComapre
This tool compares a pair of data structures [ files, directoris and archives ]
Provides concise information about the archive content using tolls defoned in the config file.

## Design

Many components of this system are heavily driven by configuration files.  This
is to allow new validation code to be added and incorporated without modifying
the driver code.

## Tools

`cgpCompare` has multiple commands, listed with `cgpCompare --help`.

### cgpCompare

Takes the input archives,files,folders and does the comparison for matching file types
based on tools defined in  `archCompare/config/*.json`
file.

Valid input types include:

* .tar - archive containing multiple files and folders to compare
* folder - any folder containing sub folders and files
* file - any file with extension configured in the `fileTypes.json` configuration file

The output is a tab separated columns containing:

* `File_a`  - compared file name  from first archive
* `File_b`  - compared file name  from second archive
* `Status`  - comparsion status [ compared, skipped ]
* `SimilarityBy` - if files are compared and found similar it will have one of the value [name , data or checksum] 
              otherwise 'differ', reason if files were skipped from comparison

Various exceptions can occur for malformed files.

## INSTALL

Installation is via `easy_install`.  Simply execute with the path to the compiled
'egg':

```bash
easy_install archCompare.egg-X.X.X-py3.6.egg

```

Installing via `pip install` .Simply execute with the path to the compiled 'whl':
```bash

pip install archCompare.X.X.X-py3-none-any.whl

```


### Package Dependancies

`easy_install` will install the relevant dependancies, listed here for convenience:

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

The release is handled by wheel:

```bash
$ source env/bin/activate # if not already
$ python setup.py bdist_wheel -d dist
# this creates an wheel archive which can be copied to a deployment location, e.g.
$ scp archCompare.X.X.X-py3-none-any.whl user@host:~/wheels
# on host
$ pip install --find-links=~/wheels archCompare
```
