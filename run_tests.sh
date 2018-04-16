#!/usr/bin/env bash
set -e
MYPACKAGE=pyCRISPRcleanR
pytest \
 --cov-report term \
 --cov-report html \
 --cov-fail-under=50 \
 --cov=$MYPACKAGE \
 -vv

#echo -e "\n#################\n# Running pycodestype:\n"
set +e
pycodestyle \
$MYPACKAGE

set +e
radon cc \
 -nc \
$MYPACKAGE
