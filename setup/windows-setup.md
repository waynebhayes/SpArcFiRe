# SpArcFiRe for Windows

## Introduction:
The SpArcFiRe repo is designed to run on Linux based computers, fortunately it can be run on Windows using WSL.

WSL (or Windows Subsystem for Linux) is a compatiablity layer that enables you to run native Linux command-line tools directly on Windows[^1]. 

In this instruction set, we will cover:
1. [Installing WSL](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#section-1-installing-wsl)
2. [Using WSL Terminal and Cloning SpArcFiRe repo](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#section-2-using-wsl-terminal-and-cloning-sparcfire-repo)
3. [Setting up SpArcFiRe for WSL](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#section-3-setting-up-sparcfire-for-wsl)
    - [Installing MATLAB Runtime](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#install-matlabs-runtime-compiler-version-r2017a)
    - [Installing Python2.7](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#install-python27-and-python-libraries-numpy-astropy-scipy-pillow)
    - [Installing GCC](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#install-gcc-a-c-compiler)
    - [Compiling Scripts](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#compile-the-c-scripts)
    - [Creating a Symbolic Link to SpArcFiRe](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#creating-a-symbolic-link-to-sparcfire)
    - [Setting up SpArcFiRe Directories](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#creating-a-symbolic-link-to-sparcfire)
    - [Downloading Sample Fits Files](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#setting-up-sparcfire-directories)
4. [Running SpArcFiRe](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#section-4-running-sparcfire)

## Section 1: Installing WSL
To install WSL, follow [the instructions here.](https://docs.microsoft.com/en-us/windows/wsl/install-win10)

**NOTE: For this tutorial, we will be using Ubuntu. You can use a different flavor of linux, but you may have to recompile the MATLAB scripts.**

## Section 2: Using WSL Terminal and Cloning SpArcFiRe repo
Now that you have installed WSL, open up the WSL terminal by clicking on your Window's Search Bar and type `Ubuntu`. The terminal will be called something like `Ubuntu 18.04.05 for Window`. Click on the icon to open the terminal. 

By default, you should be in your home directory `/home/YOURUSERNAME` (for example for me, my home directory is `/home/cora`). To double check where you are, you can type the command `pwd` which will print the working directory. If you appear to be in a different driectory, type the command `cd ~` to change your current directopry to your home directory.

Once you are in your home directory, we can clone the SpArcFiRe repo using git. You will likely need to install git, which can be done using the following commands:
```
sudo apt-get update
sudo apt-get install git -y
```
A few things to note: both these commands require sudo to run. For those not familiar, sudo is essentially the same as Administartor on Windows. The first command gets the most recent packages avalible, and the second command install git (essentially how you can use github on the command line).

After installing git, make sure you are in your home directory, and then you can clone the repo using the bellow command [^2].
```
git clone https://github.com/waynebhayes/SpArcFiRe.git
```

## Section 3: Setting up SpArcFiRe for WSL
Now that you have cloned the SpArcFiRe repo, we are ready to setup SpArcFiRe for WSL.

Now I have good news, and bad news. First the bad news, we have to all of the following: 
* Install unzip and ImageMagick
* Install MATLAB's Runtime Compiler (version R2017a)
* Install Python2.7 and Python Libraries: numpy, astropy, scipy, Pillow
* Install GCC (a C-compiler)
* Compile the C scripts
* Creating a Symbolic Link to SpArcFiRe
* Setting up SpArcFiRe Directories
* Downloading Sample Fits Files (for testing)

But the Good News: I have provided scripts to do all the steps above automatically. I have provided two versions of the instructions:
* [Full WSL Setup](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#full-setup) - which runs all setup steps (use this for a quick setup)
* [Step-Step-Setup](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#step-step-setup) - shows the inidivual set up steps (use this to pick and choose what to run/debug)

### Full Setup:
**Important Note:** By running this script you are agreeing to MATLAB's terms and conditions [^3] and copy of which can be found [here](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/licenses/MCR_license.txt).
Begin by navigating to the setup directory: `cd SpArcFiRe/setup/`
Now type the command: `./run-full-wsl-setup.bash`
If you want to install GALFIT, type the command:`./gulfit-setup.bash`
If there was an issue, take a look at the step by step instructions (you may just have to run a single step again). If everything works, congrats, now it is time to [run SpArcFiRe](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#section-4-running-sparcfire)

### Step-Step-Setup:
**Important Note:** Please make sure you are in the following directory before running any of these commands: `cd ~/SpArcFiRe/setup/setup-scripts/`

#### Install Unzip and ImageMagick
Command: `sudo apt install unzip`
Command: `sudo apt install imagemagick`
Usage: These are useful in the following installation, if you haven't install unzip, you may can't install MATLAB's Runtime Compiler successfully.

#### Install MATLAB's Runtime Compiler (version R2017a)
**Important Note:** By running this script you are agreeing to MATLAB's terms and conditions [^3] and copy of which can be found [here](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/licenses/MCR_license.txt).
Command: `install-mrc.bash`
Usage: This Install MATLAB's Runtime Compiler (version R2017a) which is used to run all the MATLAB SpArcFiRe scripts. This is done using silent install (meaning it won't need any additional input to run).

#### Install Python2.7 and Python Libraries: numpy, astropy, scipy, Pillow ,tsv
Command: `install-python2-plus-libraries.bash`
Usage: Installs Python2.7 and libraries to run SpArcFiRe's script. Note that Python2.7 has been deprecated, but as of now the code has not been moved to Python3 (so we use Python2.7).

#### Install GCC (a C-compiler)
Command: `install-gcc.bash`
Usage: Installs the C compiler gcc. This is used to compile the C scripts used in sparcfire (such as delete-commas-inside-quotes.c which is used to make tsv files)

#### Compile the C scripts
Command: `compile-scripts.bash`
Usage: Compiles the C scripts.

#### Creating a Symbolic Link to SpArcFiRe
Command: `create-symbolic-link.bash`
Usage: Create a symbolic link in your home directory called bin which links to `~/SpArcFiRe/scripts/`. This is used in many SpArcFiRe scripts to locate other code files.

#### Setting up SpArcFiRe Directories
Command:`setup-dirs.bash`
Usage: Creates the following directories:
* ~/SDSS/G.in: Where you place input files for SpArcFiRe
* ~/SDSS/G.tmp: Where temporary files generated by SpArcFiRe are placed'
* ~/SDSS/G.out: Where SpArcFiRe's output is found
* ~/SDSS/G.out/matout: Where some of SpArcFiRe's matlab output files are found

#### Downloading Sample Fits Files (for testing)
Command: `download-sample-fits.bash`
Usage: To test out SpArcFiRe you will need some input, this command downloads 3 sample fits files and places them in the input folder.

## Section 4: Running SpArcFiRe
After the [WSL Setup is complete](https://github.com/cora-schallock/SpArcFiRe/blob/master/setup/windows-setup.md#section-3-setting-up-sparcfire-for-wsl), we are now ready to run SpArcFiRe.

First navigate to the SpArcFiRe directory and check your setup is valid:
```
cd ~/SpArcFiRe
source ./setup.bash ~/SpArcFiRe
```
If you want to run SpArcFiRe with the defult command line arguments, you can run SpArcFiRe with the command: `~/scripts/wschallo/SpArcFiRe-run.sh`

This command is equiallent to running:
"$HOME"/SpArcFiRe/scripts/SpArcFiRe -convert-FITS $HOME/SDSS/G.in $HOME/SDSS/G.tmp $HOME/SDSS/G.out  -generateFitQuality 0 -writeBulgeMask 1

If you want to run SpArcFiRe with different command line arguments, please see our [README.md](https://github.com/waynebhayes/SpArcFiRe#readme)).

## Section 5: Install package needed for GALFIT
All the effort above is to make sure you can run the basic function of SpArcFiRe, and if you want to run GALFIT module, these following packages are needed for python3

Install numpy,if you haven't, by using Command `python3 -m pip instal1 numpy`
Install scipy,if you haven't, by using Command `python3 -m pip instal1 scipy`
Install astropyif you haven't, by using Command `python3 -m pip install astropy`
Install pandas if you haven't, by using Command `python3 -m pip install pandas`
Install lpython if you haven't, by using Command `python3 -m pip instal1 Ipython`
Install skimage if you haven't, by using Command `python3 -m pip install scikit-image`
Install matplotlib if you haven't, by using Command `python3 -m pip instal1 matplotlib`
Install pickle if you haven't by using Command `python3 -m pip install pickle-mixin`
Install joblib if you haven't, by using Command `python3 -m pip install joblib`

By the way, we also need to install GALFIT, if you want to learn more about it, click [here](https://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html) to get more information. Here are instructions about how to install galfit.

1) Download [GALFIT for Ubuntu](https://users.obs.carnegiescience.edu/peng/work/galfit/galfit3-debian64.tar.gz) and remove the package to `~/.local/bin`
2) unzip the package by command `tar -zxvf galfit3-debian64.tar.gz`(change the file name if necessary)
3) delete the package by `rm galfit3-debian64.tar.gz`
4) use command `ls` to make sure galfit it in the list

## Section 6: Test SpArcFiRe
Now you finall have installed everything needed for SpArcFiRe, try to run the test to check your efforts!
```
cd ~/SpArcFiRe
./regression-test-all.sh
```

## Section 7: Debugging Issues
Please refer to the debugging list in the [windows-debug.md](https://github.com/waynebhayes/SpArcFiRe/blob/master/setup/windows-debug.md) file.

## References:
[^1]: [WSL](https://docs.microsoft.com/en-us/windows/wsl/install)
[^2]: [Cloning a Github Repo](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
[^3]: [MATLAB](https://www.mathworks.com/help/compiler_sdk/dotnet/install-the-matlab-runtime.html)
