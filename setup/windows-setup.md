# SpArcFiRe for Windows

## Introduction:
The SpArcFiRe repo is designed to run on Linux based computers, fortunately it can be run on Windows using WSL.

WSL (or Windows Subsystem for Linux) is a compatiablity layer that enables you to run native Linux command-line tools directly on Windows[^1]. 

In this instruction set, we will cover:
1. Installing WSL
2. Using WSL Terminal and Cloning SpArcFiRe repo
3. Setting up SpArcFiRe for WSL
    - Installing MATLAB Runtime
    - Installing Python2.7
    - Installing GCC
    - Compiling Scripts
    - Creating a Symbolic Link to SpArcFiRe
    - Setting up SpArcFiRe Directories
    - Downloading Sample Fits Files
4. Running SpArcFiRe

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

After installing git, you can clone the repo [^2].

## References:
[^1]: [WSL](https://docs.microsoft.com/en-us/windows/wsl/install)
[^2]: [Cloning a Github Repo](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
