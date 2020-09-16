#!/bin/sh
echo "\n\n---Gettting latest package lists---\n\n"
apt-get update
apt-get install unzip

echo "\n\n---Installing GIT---\n\n"
apt-get install git -y

echo "\n\n---Installing GCC---\n\n"
apt install gcc -y


echo "\n\n---Installing Python2.7, numpy, astropy, scipy, and Pillow---\n\n"
apt install python2.7 python-pip -y
pip2 install numpy
pip2 install astropy
pip2 install scipy
pip2 install Pillow

# create a temporary directory to download and unzip matlab installs
mkdir /matlab_downloads/

echo "\n\n---Downloading MATLAB Compiler 2017a---\n\n"
wget -P /matlab_downloads/ https://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_glnxa64_installer.zip
unzip /matlab_downloads/MCR_R2017a_glnxa64_installer.zip -d /matlab_downloads/
echo "\n\nInstalling MATLAB Compiler 2017a."
/matlab_downloads/install -mode silent -agreeToLicense yes -destinationFolder /pkg/matlab

# cleanup all temp directories
rm -r /matlab_downloads/*
rmdir /matlab_downloads

# rename directory to the one SpArcFiRe expects
mv /pkg/matlab/v92 /pkg/matlab/R2017a

echo "\n\n---Cloning SpArcFiRe---\n\n"
sudo -u username git clone https://github.com/wschallo/SpArcFiRe.git


echo "\n\n---Setting Up SpArcFiRe---\n\n"
sudo -u username SpArcFiRe/scripts/wschallo/setup.sh

echo "Setup Complete, to run SpArcFiRe on the test images run ~/bin/wschallo/SpArcFiRe-run.sh"
echo "For an in-depth tutorial on how to use SpArcFiRe read the README file or look at it on the GitHub repo."
