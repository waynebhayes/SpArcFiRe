#Description:
# SpArcFiRe with/without Galfit comparison
#Usage:
# - Run on ics.uci.edu servers to access sparcfire

# Activate python virtualenv
source /home/dlcheng/sparcfire/scripts/SpArcFiRe-pyvenv/bin/activate

#Run SpArcFiRe on images with/without Galfit
#echo "Running SpArcFiRe without Galfit"
echo "Running SpArcFiRe with Galfit"
/home/dlcheng/bin/SpArcFiRe -convert-FITS /home/dlcheng/SDSS/G.in /home/dlcheng/SDSS/G.tmp /home/dlcheng/SDSS/G.out -generateFitQuality 1 -writeBulgeMask 1

# Generate images for html page
./generate_html.sh

# Deactivate python virtualenv
deactivate
