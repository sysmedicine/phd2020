# Workshop Preparation

## Introduction
In order to make sure that we have a smooth and on-time workshop sessions, please make sure that you do all the steps below before the course

**We will have a Lab Prep Session on the first day, but please do this before the course. The purpose of the session is to fix or answers any problems or questions of the setup**

**Any problem that you find prior to the course, please email muhammad.arif [at] scilifelab.se**

## Preparation Steps
We will run most of our softwares inside a conda environment. You can either download the conda environment or download the pre-installed virtual machine (with virtual box)

### Option 1: Download Conda Environment
This step is the lightweight step with relatively less usage of your space. But due to the complications with different OS, we would recommend this if you are familiar with conda or unix, with Linux/MacOSX computers. Otherwise, option 2 is the best options
1. Install Miniconda in your computer. To install miniconda:
> #MacOSX  
> curl -o Miniconda3-latest-MacOSX-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh  
> sh Miniconda3-latest-MacOSX-x86_64.sh  

> #Linux  
> wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  
> sh Miniconda3-latest-Linux-x86_64.sh  

> #Windows  
> Our suggestion is to use Option 2.

2. Download conda environment yaml file
> #MacOSX  
> curl -o env.yaml https://raw.githubusercontent.com/sysmedicine/phd2020/master/LabPrep/env.yaml

> #Linux  
> wget https://raw.githubusercontent.com/sysmedicine/phd2020/master/LabPrep/env.yaml

3. Create a new conda environment with the file
> conda env create -p phd2020 -f env.yaml


