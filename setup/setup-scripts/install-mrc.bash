#!/bin/bash

# IMPORTANT NOTE!! : By running this script, you will agree to mrc's license
# NOTE: for more info about setting up mrc using silent mode, please see: https://www.mathworks.com/help/compiler_sdk/dotnet/install-the-matlab-runtime.html

#1) Create tmp dir for matlab compiler install:
mkdir matlab_tmp

#2) download mrc r2017a to tmp dir:
wget https://ssd.mathworks.com/supportfiles/downloads/R2017a/deployment_files/R2017a/installers/glnxa64/MCR_R2017a_glnxa64_installer.zip -P matlab_tmp

#3) unzip mrc:
unzip matlab_tmp/MCR_R2017a_glnxa64_installer.zip -d matlab_tmp

#4) silent setup:
## source about silent setup params: https://www.mathworks.com/help/compiler_sdk/dotnet/install-the-matlab-runtime.html
sudo ./matlab_tmp/install -mode silent -agreeToLicense yes -destinationFolder /pkg/matlab

#5) move matlab:
sudo mv /pkg/matlab/v92 /pkg/matlab/R2017a

#6) Delete tmp dir:
rm -rf matlab_tmp
