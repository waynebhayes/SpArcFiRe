#!/bin/bash

#0) Get latests packages:
sudo apt-get update

#1) Install MRC:
sudo ./setup-scripts/install-mcr.bash

#2) Install Python 2.7 plus Libraries:
sudo ./setup-scripts/install-python2-plus-libraries.bash

#3) Install gcc:
sudo ./setup-scripts/install-gcc.bash

#4) Create symbolic link: "$HOME"/bin -> SpArcFiRe/scripts
sudo ./setup-scripts/create-symbolic-link.bash

#5) Create dir structure:
./setup-scripts/setup-dirs.bash

#6) Download sample fits:
./setup-scripts/download-sample-fits.bash
