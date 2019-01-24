#!/bin/bash

#create output directories
echo "$HOME"
mkdir "$HOME"/SDSS
mkdir "$HOME"/SDSS/G.in
mkdir "$HOME"/SDSS/G.out
mkdir "$HOME"/SDSS/G.tmp

#download sample fits files to test SpArcFiRe on
wget -O "$HOME"/SDSS/G.in/1237648702966988820_g.fits http://sparcfire.ics.uci.edu/sampleFitsFiles/1237648702966988820_g.fits
wget -O "$HOME"/SDSS/G.in/1237648702967054406_g.fits http://sparcfire.ics.uci.edu/sampleFitsFiles/1237648702967054406_g.fits
wget -O "$HOME"/SDSS/G.in/1237648702967251093_g.fits http://sparcfire.ics.uci.edu/sampleFitsFiles/1237648702967251093_g.fits

#make symbolic link to access SpArcFiRe
ln -s "$HOME"/SpArcFiRe/scripts/ "$HOME"/bin

#edit Path variables
sed -i 's#/home/dlcheng/sparcfire/scripts#"$HOME"/bin#g' your_file
sed -i 's/dlcheng/wschallo/g' "$HOME"/bin/SpArcFiRe
