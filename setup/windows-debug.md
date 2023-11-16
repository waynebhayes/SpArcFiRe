
# SpArcFiRe for Windows - Debugging

## Introduction:
This document contains a list of current issues with running SpArcFiRe on WSL and information about the steps taken to fix the issues.


## Section 1: Current Issues
#### Issue #1: Matlab mrc r2017a shipping with outdated .so files
If you encounter an error along the lines of:
 `/home/cora/SpArcFiRe/scripts/ArcServer/findClusterArcsServer.r140: symbol lookup error: /lib/x86_64-linux-gnu/libfontconfig.so.1: undefined symbol: FT_Done_MM_Var` 
The libfreetype dynamic library file is out of date.
Solution: [Set aside libfreetype.so.6](https://www.mathworks.com/matlabcentral/answers/364551-why-is-matlab-unable-to-run-the-matlabwindow-application-on-linux#answer_288902)

Step-by-Step Guide:
1) Use the command `locate libfreetype` to locate the libfreetype dynamic library file 
*Note: if command locate not found, install locate with `sudo apt install plocate`*
You are looking for the path of a libfreetype.so.6 inside your matlab bin directory. For me it was located here: `/pkg/matlab/R2017a/bin/glnxa64/libfreetype.so.6`
2) cd to the directory of the file found in step 1
For example: `/pkg/matlab/R2017a/bin/glnxa64`
3) Exclude the file by creating an excluded directory and moving the old dynamic library file to the excluded directory:
`sudo mkdir exclude`
`sudo mv libfreetype.so* exclude/`
4) Try Rerunning SpArcFiRe

Sources:
[Matlab Answers](https://www.mathworks.com/matlabcentral/answers/364551-why-is-matlab-unable-to-run-the-matlabwindow-application-on-linux#answer_288902)
[Reddit post about Matlab breaking](https://www.reddit.com/r/archlinux/comments/tkas9q/matlab_stopped_working_after_system_update/)

#### Issue #2: libXt.so.6: cannot open shared object file: No such file or directory
If you encounter an error along the lines of: `libXt.so.6: cannot open shared object file: No such file or directory` or `libstdc++.so.6: cannot open shared object file: No such file or directory`

You will need to fix an issue with the `libXt.so.6`/`libstdc++.so.6` files included in Matlab.
Solution: [Follow this Windows Matlab Setup](https://www.mathworks.com/matlabcentral/answers/308911-can-i-install-matlab-in-bash-on-ubuntu-on-windows#answer_333214)
Step-by-Step Guide:
1) Update basic library:
```
sudo apt-get install libstdc++6
sudo add-apt-repository ppa:ubuntu-toolchain-r/test 
sudo apt-get update
sudo apt-get upgrade
sudo apt-get dist-upgrade
```

2) If you are using [opengl](https://en.wikipedia.org/wiki/OpenGL) install a display library like mesa librar:
`sudo apt-get install mesa-utils`

3) Install jdk:
`sudo apt-get install default-jdk`

4) Install execstack:
`sudo apt-get install execstack`

5) Use the command `locate libfreetype` to locate the libfreetype dynamic library file 
*Note: if command locate not found, install locate with `sudo apt install plocate`*
You are looking for the path of a libfreetype.so.6 inside your matlab bin directory. For me it was located here: `/pkg/matlab/R2017a/bin/glnxa64/libfreetype.so.6`

6) Change directories: `cd /pkg/matlab/R2017a/bin/glnxa64`

7) Call execstack:
```
execstack -c libmwblas.so
execstack -c libmwlapack.so
```
8) Rerun SpArcFiRe

Sources:
[Matlab Answers](https://www.mathworks.com/matlabcentral/answers/308911-can-i-install-matlab-in-bash-on-ubuntu-on-windows#answer_333214)

#### Issue #3: Abnormal termination: Segmentation violation

*Note: This is similar to the libstdc++.so.6 in issue #2 (that solution may work to resolve this issue)*

If you encounter an error along the lines of: 
```
Abnormal termination:

Segmentation violation
Register State (from fault):

RAX = 0000000000000000 RBX = 00007fdb91e76808
...
```
This maybe an issue with installing matlab on a newer version of Ubuntu (or a new version of Ubuntu on WSL).
Solution: [Follow this guide](https://www.mathworks.com/matlabcentral/answers/275176-matlab-crashes-on-startup-segmentation-violation#answer_231458)

Step-by-Step Guide:
1) Locate the c++ .so file with command `locate libstdc++.so.6`
You should see a path that looks like: `/pkg/matlab/R2017a/sys/os/glnxa64/libstdc++.so.6` and/or `/pkg/matlab/R2017a/sys/os/glnxa64/libstdc++.so.6.0.20`

2) Navigate to the directory from step 1. For example: `cd /pkg/matlab/R2017a/sys/os/glnxa64`

3) Rename the libstdc++.so.6 files. For example: `sudo mv libstdc++.so.6 libstdc++.so.6.old` or `sudo mv libstdc++.so.6.0.20 libstdc++.so.6.0.20.old`
4) Rerun SpArcFiRe

Sources:
[Matlab Answers](https://www.mathworks.com/matlabcentral/answers/275176-matlab-crashes-on-startup-segmentation-violation#answer_231458)


#### Issue #4: libncurses.so.5: cannot open shared object file: No such file
Note: This is similar to the libstdc++.so.6 in issue #2 and #3
If you encounter an error along the lines of: `libncurses.so.5: cannot open shared object file: No such file`, follow these commands may help
```
sudo apt install apt-file
sudo apt-file update
sudo apt-file find libncurses.so.5
sudo apt install libncurses5
sudo apt install libncurses*
```

Sollution:[Follow this guide](https://blog.csdn.net/qq_35078688/article/details/125326873)

#### Issue #4: Warning: header column name "XXX" does not appear in all input files
if you face the warning like this or like below, it may result in lacking tsv module for python2
```
Warning:
****************
****************
**************** ABOVE ERRORS (60 header, 0 diff) MAY BE FATAL; SUPPRESSING FURTHER WARNINGS
****************
****************
```
Solve this by command `python2 -m install tsv`