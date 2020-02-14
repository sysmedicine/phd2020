##  Preparation Step for GEM

To begin our lab, we need to install all the needed packages to prepare for the GEM based analysis. COBRA and RAVEN toolboxes are the two main toolboxes for GEM based analysis, and in this lab, we will be mainly focusing on RAVEN toolbox. If you are interested in COBRA toolbox, you could go to check [the link](https://opencobra.github.io/cobratoolbox/stable/index.html).

For installation of RAVEN, you should first go to [this link](https://github.com/SysBioChalmers/RAVEN) and download the package as a zip file, and unzip it to the path you defined. After that, you should include this path in MATLAB to enable the usage of RAVEN functions.

The RAVEN toolbox is dependent on the Systems Biology Markup Language (SBML) library for importing/exporting GEMs in SBML format. The libSBML MATLAB API 5.17 is recommended, and you can go [this link](https://sourceforge.net/projects/sbml/files/libsbml/5.17.0/stable/MATLAB%20interface/) to download the library API for MATLAB. Again, you will need to unzip the downloaded package to a path you select, and include the path in MATLAB.

In addition, linear programming solver is needed for the model simulation for RAVEN toolbox. In this lab, we will use MOSEK for the modeling. You can download it (version 7.1.0.63) through (the link)[https://www.mosek.com/downloads/7.1.0.63/], and install it follow the instruction (Note that if you are using a computer in the lab, it should have been already installed). You will also need to apply for an academic license from [the website](https://www.mosek.com/products/academic-licenses/), and you will get instant automatic response if you apply with your academic email. After installation and license the MOSEK solver following the instruction in the email, you will need to link the solver to MATLAB by including its path.