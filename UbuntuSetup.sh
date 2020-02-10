#!/bin/sh
echo "\n\n---Gettting latest package lists---\n\n"
sudo apt-get update

echo "\n\n---Installing GIT---\n\n"
sudo apt-get install git -y

echo "\n\n---Installing GCC---\n\n"
sudo apt install gcc -y


echo "\n\n---Installing Python2.7, numpy, astropy, scipy, and Pillow---\n\n"
sudo apt install python2.7 python-pip -y
pip2 install numpy
pip2 install astropy
pip2 install scipy
pip2 install Pillow

# create a temporary directory to download and unzip matlab installs
sudo mkdir /matlab_downloads/

echo "\n\n---Downloading MATLAB Compiler 2014a---\n\n"
sudo wget -P /matlab_downloads/ https://www.mathworks.com/supportfiles/downloads/R2014a/deployment_files/R2014a/installers/glnxa64/MCR_R2014a_glnxa64_installer.zip
sudo unzip /matlab_downloads/MCR_R2014a_glnxa64_installer.zip -d /matlab_downloads/
echo "\n\nAbout to install MATLAB Compiler 2014a,"
echo "IMPORTANT: When prompted set the installation folder to /pkg/matlab"
echo "			 This must be done correctly for SpArcFiRe to work."
echo "Press ENTER when ready to install."
read -p "" key
sudo /matlab_downloads/install
# clear out the temp directory for the next install
sudo rm -r /matlab_downloads/*

echo "\n\n---Downloading MATLAB Compiler 2017a---\n\n"
sudo wget -P /matlab_downloads/ https://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_glnxa64_installer.zip
sudo unzip /matlab_downloads/MCR_R2017a_glnxa64_installer.zip -d /matlab_downloads/
echo "\n\nAbout to install MATLAB Compiler 2017a."
echo "IMPORTANT: When prompted set the installation folder to /pkg/matlab"
echo "			 This must be done correctly for SpArcFiRe to work."
echo "Press ENTER when ready to install."
read -p "" key
sudo /matlab_downloads/install

# cleanup all temp directories
sudo rm -r /matlab_downloads/*
sudo rmdir /matlab_downloads

echo "\n\n---Cloning SpArcFiRe---\n\n"
git clone https://github.com/wschallo/SpArcFiRe.git


echo "\n\n---Setting Up SpArcFiRe---\n\n"
SpArcFiRe/scripts/wschallo/setup.sh

echo "Setup Complete, to run SpArcFiRe on the test images run ~/bin/wschallo/SpArcFiRe-run.sh"
echo "For a more in-depth tutorial on how to use SpArcFiRe, visit the github page under intructions."


# FILES THAT HAVE TO BE CHANGED FOR THIS TO WORK
# /bin/ArcPNGcreator/convertFromFits.0.sh - set MCRROOT=/pkg/matlab/v83, set R=121
# /bin/ArcsServer/run_findClusterArcsServer.sh - set MCRROOT=/pkg/matlab/v92
# /bin/arc2csv.sh - set MCRROOT=/pkg/matlab/v92
# /bin/writeGxyParamsToCsv.sh - set MCRROOT=/pkg/matlab/v92
# /SpArcFiRe/SpArcFiRe - replace all "python -s" with "python2.7", line 205 == should be =