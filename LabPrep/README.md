# Workshop Preparation

## Introduction
In order to make sure that we have a smooth and on-time workshop sessions, please make sure that you do all the steps below before the course

**Installation might take a while, so please do it ASAP! If you have any permission problems, please contact your IT team, since we will not be able to help you on this. We will have a Lab Prep Session on the first day, but please do this before the course. The purpose of the session is to fix or answers any problems or questions of the setup**

**Any problem that you find prior to the course, please email muhammad.arif [at] scilifelab [dot] se**

## Preparation Steps
We will run most of our softwares inside a conda environment. You can either download the conda environment or download the pre-installed virtual machine (with virtual box). Please pick one of out of three options below.

### Option 1: Using Virtual Machine with GUI (Complexity: Lowest; Download Size: Highest)
Download Size: 8 Gb

We require around 25 Gb of empty space in harddrive and mininum 4 Gb RAM. Higher RAM is definitely better.

1. Download Virtual Box: https://www.virtualbox.org/wiki/Downloads
2. Download the pre-installed Ubuntu virtual machine with GUI on [this link](https://kth.box.com/s/4u1m5tjmv3eeh3q4v7vr8innle13cxcx)
3. Click on the option "File" --> "Host-network manager" in the toolbar and click "Create" until you see a "vboxnet0" in the table.
4. Import the virtual machine to Virtual Box
5. When needed, use sysmedicine as both username and password to log in to the virtual machine.

With virtual machine, there are several ways to access the terminal:
1. Directly from virtual box, there's an installed GUI.
2. "Remotely" via SSH --> (can be done via terminal in Linux/MacOSX, or MobaXterm in Windows)

Inside the virtual machine, all the necessary software and GUI has been installed. In order to access it remotely, check the IP address of the virtual machine by "ifconfig | grep inet" and look for the IP address 192.xxx.xxx.xxx

1. To access Rstudio: http://192.xxx.xxx.xxx:8787
2. (After typing "jupyterhub" in terminal) To access Jupyter Notebook: http://192.xxx.xxx.xxx:8000

**EncyclopeDIA must be run from inside the virtual machine with GUI**

**OBS: Don't update the OS or the conda packages unless instructed by the teachers.**

### Option 2: Using Virtual Machine no GUI (Complexity: Intermediate; Download Size: Intermediate)
Download Size: 4.4 Gb

We require around 15 Gb of empty space in harddrive and mininum 4 Gb RAM. Higher RAM is definitely better.

1. Download Virtual Box: https://www.virtualbox.org/wiki/Downloads
2. Download the pre-installed Ubuntu virtual machine without GUI on [this link](https://kth.box.com/s/isnho6tr65omzvxcdnrqp85gyuyv84zc)
3. Click on the option "File" --> "Host-network manager" in the toolbar and click "Create" until you see a "vboxnet0" in the table.
4. Import the virtual machine to Virtual Box
5. When needed, use sysmedicine as both username and password to log in to the virtual machine.
6. Since there's no GUI in the server, EncyclopeDIA has to be installed manually in your own computer. Please follow [this link](https://bitbucket.org/searleb/encyclopedia/wiki/Home) to install EncyclopeDIA

The virtual machine can be acccessed "Remotely" via SSH --> (can be done via terminal in Linux/MacOSX, or MobaXterm in Windows)

Inside the virtual machine, all the necessary softwares have been installed. In order to access it remotely, check the IP address of the virtual machine by "ifconfig | grep inet" in terminal inside the VM and look for the IP address 192.xxx.xxx.xxx

1. To access Rstudio: http://192.xxx.xxx.xxx:8787
2. (After typing "jupyterhub" in terminal) To access Jupyter Notebook: http://192.xxx.xxx.xxx:8000

**OBS: Don't update the OS or the conda packages unless instructed by the teachers.**

### Option 3: Download Conda Environment (Complexity: Highest;Download Size: Lowest)
This step is the lightweight step with relatively less usage of your space. But due to the complications with different OS, we would recommend this if you are familiar with conda or unix, with Linux/MacOSX computers. Otherwise, option 1 is the best options

The steps were validated to work properly in the OS below:
1. MacOSX Catalina
2. Ubuntu 18.04 LTS

**OBS: Based on past experience, there were too many problems between conda and Windows, especially the newer Windows versions. If you  want to use conda in Windows, do it at your own risk. We will not be supporting any questions about how to use the conda environment inside windows. Our answer will only be: use the provided virtual machine in Option 2**

**If you would like to use Conda Environment on Windows 10, you can refer to [this link](https://nbisweden.github.io/workshop-scRNAseq/conda_instructions.html), and use our environment yaml file instead of theirs**

1. Install Miniconda in your computer. To install miniconda:
> #MacOSX  
> curl -o Miniconda3-latest-MacOSX-x86_64.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh  
> sh Miniconda3-latest-MacOSX-x86_64.sh  

> #Linux  
> wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  
> sh Miniconda3-latest-Linux-x86_64.sh  

> #Windows  
> Our suggestion is to use Option 2. 

2. Restart your terminal

3. Download conda environment yaml file
> #MacOSX  
> curl -o env.yaml https://raw.githubusercontent.com/sysmedicine/phd2020/master/LabPrep/env.yaml

> #Linux  
> wget https://raw.githubusercontent.com/sysmedicine/phd2020/master/LabPrep/env.yaml

4. Create a new conda environment with the file
> conda env create -f env.yaml

**OBS: In some computers, the packages installed under #Metagenomics can cause a problem. Our advice is if you have a persistent error, try to comment/remove all the packages under #Metagenomics in the yaml file.**

5. Install R and Rstudio (you can use conda or standalone). For Ubuntu, please follow [this link](https://www.digitalocean.com/community/tutorials/how-to-install-r-on-ubuntu-18-04) to install the latest R. Make sure that you have the latest version (3.6)

6. (Ubuntu Only) Install with apt-get several libraries
> sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev

7. Download R script for the installation
> #MacOSX  
> curl -o env.yaml https://raw.githubusercontent.com/sysmedicine/phd2020/master/LabPrep/R_install_pkgs.R

> #Linux  
> wget https://raw.githubusercontent.com/sysmedicine/phd2020/master/LabPrep/R_install_pkgs.R

8. Run the script to install R packages. You can do it via Rstudio.

9. Install MinPath by following [this link](https://omics.informatics.indiana.edu/MinPath/)

10. [Install SHORTBRED](https://huttenhower.sph.harvard.edu/shortbred). Don't forget to get the license for USEARCH. Blast, Muscle, CD-HIT and Biopython requirement has already been installed in the previous step.

## Conda Environment

For all options, some packages are installed under conda environment. By default, the name of the environment is "phd2020". To activate the environment, type:
> conda activate phd2020

Until you see something like this in your terminal:
> **(phd2020)** sysmedicine@sysmedicine-vm $ 

## Problems During Installation

If you follow the steps above, hopefully no problem will arise. In case of error, please first consult Google. If there's no solution to your problem, email muhammad.arif [at] scilifelab [dot] se the error message. If problem persists, please consider using the virtual machine. 
