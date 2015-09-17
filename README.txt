====================================
Iterative Feature Removal
Source code written by:
Stephen O'Hara, Colorado State University
Copyright 2012-2013
====================================

This is open-source software provided to allow
other researchers to work with Iterative Feature Removal (IFR),
and to reproduce the results described in our paper:

Stephen O'Hara, Kun Wang, Richard A. Slayden, Alan R. Schenkel, Greg Huber, 
Corey S. O'Hern, Mark D. Shattuck and Michael Kirby
"Iterative Feature Removal Yields Highly Discriminative Pathways,"
BMC Genomics, under review, 2013

############################
Installation
############################
The code provided in this package is dependent upon a few open source libraries
beyond those that are core to scientific programming in python. You will need
python 2.7, with numpy/scipy and matplotlib. You will also need scikit-learn
for python 2.7. Google it, and you will find information on how to install it.

In addition, to produce LaTeX format output for the tables, we use a package
called asciitable. Again, if you search for it, you can find it. This README
is not a primer on how to install and maintain open source software.

If you use a macintosh computer, all required packages are easily installed via
macports.

############################
# GETTING STARTED          #
############################
There are two mostly-equivalent tutorials in the \doc folder. They
are 'tutorial.txt', which is plain-text explanation of how to get
started, how to run scripts to reproduce the figures and tables from
our paper, and how one can use IFR on one's own data.

The second tutorial is an interactive version using iPython Notebook.
If you have an up-to-date version of iPython installed with the notebook
functionality, you can work with this tutorial as follows.

1. Go to a command prompt and navigate to the directory (...\IFR_Dist\doc)
where the notebook file is stored. "IFR Tutorial.ipynb"

2. Launch the notebook server from the command prompt
> ipython notebook --pylab=inline

3. It should launch a web browser page, and you should see the notebook
listed in the menu. Open that notebook and you will see code and explanatory
text interwoven.

4. You can execute the contents of each cell, one-by-one, and you can try
altering the code to see what happens.

############################
# DATA ATTRIBUTION         #
############################
We employ four data sets in the analysis presented in the above paper,
and as implemented in the script file: scripts/bmc_genomics_paper.py.

Three are from the Kent Ridge Biomedical Data Repository, and the scripts
are designed to download and install them the first time they are needed.
The data repository is at: http://levis.tongji.edu.cn/gzli/data/mirror-kentridge.html
which includes citations to the original sources (as we do in our paper).

One is the Influenza data from Duke. The source of this data is available
as MATLAB files from: http://people.ee.duke.edu/lcarin/reproduce.html
We needed to re-distribute the relevant data that has been converted to
a format readable by the python/scipy scripts. The converted data
files are available in the directory: data/Duke_Data_For_Python. They
are still matlab files, but saved in an older format and restructured for
convenience. Newer matlab files cannot be read by python or other
open-source languages.

############################
# LICENSE INFORMATION      #
############################
This software is licensed under a Creative Commons license,
Attribution-NonCommercial-ShareAlike 3.0 Unported. The
details of the license are provided in the license.html file that can
be found in the same directory as this README file.

This license does not permit commercial use of this software.
Inquiries for a license allowing commercial usage should be directed to
Professor Michael Kirby at Colorado State University.