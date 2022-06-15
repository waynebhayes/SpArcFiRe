#!/bin/bash

#0) Get latests packages:
sudo apt-get update

#1) Install MRC:
sudo ./setup-scripts/install-mrc.bash

#2) Install Python 2.7 plus Libraries:
sudo ./setup-scripts/install-python2-plus-libraries.bash

#3) Install gcc:
sudo ./setup-scripts/install-gcc.bash

#4) Compile c scripts:
./setup-scripts/compile-scripts.bash

#5) Create symbolic link: "$HOME"/bin -> SpArcFiRe/scripts
sudo ./setup-scripts/create-symbolic-link.bash

#6) Create dir structure:
./setup-scripts/setup-dirs.bash

#7) Download sample fits:
./setup-scripts/download-sample-fits.bash
