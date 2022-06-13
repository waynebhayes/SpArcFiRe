#!/bin/bash

# IMPORTANT NOTE!! : By running this script, you will agree to mrc's license
# NOTE: for more info about setting up mrc using silent mode, please see: https://www.mathworks.com/help/compiler_sdk/dotnet/install-the-matlab-runtime.html

#1) Create a temporary directory to download and unzip matlab installs
sudo mkdir /matlab_downloads/

#2) Download mrc2017a to tmp dir
echo "\n\n---Downloading MATLAB Compiler 2017a---\n\n"
sudo wget -P /matlab_downloads/ https://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_glnxa64_installer.zip
sudo unzip /matlab_downloads/MCR_R2017a_glnxa64_installer.zip -d /matlab_downloads/

#3) Silent Install:
## source about silent setup params: https://www.mathworks.com/help/compiler_sdk/dotnet/install-the-matlab-runtime.html
sudo ./matlab_downloads/install -mode silent -agreeToLicense yes -destinationFolder /pkg/matlab
echo "Installed to: /pkg/matlab"

#4) Cleanup all temp directories
sudo rm -r /matlab_downloads/*
sudo rmdir /matlab_downloads

#5) Rename directory to the one SpArcFiRe expects
sudo mv /pkg/matlab/v92 /pkg/matlab/R2017a
