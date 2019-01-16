# CHANGES
## 2.0.8
* added version tag to results module, avoids failing tests

## 2.0.7
* Signture input can be directory or archive file, modified to handle docstore json input
* Signature names now automatically constructed based on file names
## 2.0.6
* Added explicit file handler close to avoid empty files in the results archive

## 2.0.5
* Updated version tag in setup.py
## 2.0.4
 * Moved tar generation step at the end to include result.html summary file in the archive
## 2.0.3
* Added utf-8 default encoding to read json, fixes text encoding issue in docker build
## 2.0.2
* renamed test file to avoid ignored by gitignore file
## 2.0.1
* remove exit statement in BAGEL
* added additional tests for bagel
## 2.0.0
* Added parallel version of BAGEL
* Added MAGeCK
* Html results summary and links to download data
* GPL3 LICENSE added
* Additional plotly figures

## 1.1.2
 * updated appropriate types [int, str etc.,] for commandline inputs

## 1.1.1
 * Added segmentation module in setup.py
 * added requirements.txt file
 * Corrected broken links in README

## 1.1.0
 * Additional output formats for normalised and corrected counts
 * segmentation is now optional

## 1.0.1
 * Updated Readme

## 1.0.0
 * First working python release of CRISPRclearR : https://github.com/francescojm/CRISPRcleanR
 * Supports parallel data analysis
 * produce interactive plotly graphs for QC analysis
