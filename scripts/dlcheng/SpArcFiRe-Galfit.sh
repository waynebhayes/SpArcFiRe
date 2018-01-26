#Description:
# SpArcFiRe with/without Galfit comparison
#Usage:
# - Run on ics.uci.edu servers to access sparcfire

# Copy images to input directory
./copy_SDSS_fits_from_file.sh 0 100

# Activate python virtualenv
source /home/dlcheng/sparcfire/scripts/SpArcFiRe-pyvenv/bin/activate

#Run SpArcFiRe on images with/without Galfit
#echo "Running SpArcFiRe without Galfit"
echo "Running SpArcFiRe with Galfit"
#Bad results - not many fits /home/dlcheng/bin/SpArcFiRe -convert-FITS /home/dlcheng/SDSS/G.in /home/dlcheng/SDSS/G.tmp /home/dlcheng/SDSS/G.out -generateFitQuality 1 -writeBulgeMask 1 -medFiltRad 8 -stopThres 0.10 -clusSizeCutoff 200
/home/dlcheng/sparcfire/scripts/SpArcFiRe -convert-FITS /home/dlcheng/SDSS/G.in /home/dlcheng/SDSS/G.tmp /home/dlcheng/SDSS/G.out -generateFitQuality 1 -writeBulgeMask 1 -medFiltRad 2 -stopThres 0.10

# Generate images for html page
./generate_html.sh

# Deactivate python virtualenv
deactivate
