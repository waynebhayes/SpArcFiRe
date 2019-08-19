#!/bin/bash

#create output directories
mkdir "$HOME"/SDSS
mkdir "$HOME"/SDSS/G.in
mkdir "$HOME"/SDSS/G.out
mkdir "$HOME"/SDSS/G.out/matout
mkdir "$HOME"/SDSS/G.tmp

#download sample fits files to test SpArcFiRe on
wget -O "$HOME"/SDSS/G.in/1237648702966988820_g.fits http://sparcfire.ics.uci.edu/sampleFitsFiles/1237648702966988820_g.fits
wget -O "$HOME"/SDSS/G.in/1237648702967054406_g.fits http://sparcfire.ics.uci.edu/sampleFitsFiles/1237648702967054406_g.fits
wget -O "$HOME"/SDSS/G.in/1237648702967251093_g.fits http://sparcfire.ics.uci.edu/sampleFitsFiles/1237648702967251093_g.fits

#make symbolic link to access SpArcFiRe
ln -s "$HOME"/SpArcFiRe/scripts/ "$HOME"/bin

#edit Path variables
#sed -i "s#/home/dlcheng/sparcfire/scripts#$HOME/bin#g" "$HOME"/bin/SpArcFiRe-pyvenv/bin/activate
sed -i "s#/home/dlcheng#$HOME#g" "$HOME"/bin/SpArcFiRe-pyvenv/bin/activate
 #  ^above line changes line: VIRTUAL_ENV="/home/dlcheng/bin/SpArcFiRe-pyvenv" in file /bin/SpArcFiRe-pyvenv/bin/activate
NAME="$(echo "$HOME" | cut -d'/' -f3)" #assume HOME in form of: /home/wschallo 
sed -i "s/dlcheng/$NAME/g" "$HOME"/bin/SpArcFiRe

#fix permission issues:
chmod +x "$HOME"/bin/ArcServer/run_findClusterArcsServer.sh
chmod +x "$HOME"/bin/wschallo/SpArcFiRe-run.sh
