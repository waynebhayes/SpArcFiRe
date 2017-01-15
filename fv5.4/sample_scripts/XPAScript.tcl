########################################################################
#
# Sample Tcl script for fv 2.6
#
# This script will open 2 sample files (available online and distributed
# with the fv executables) and displays their contents in several forms...
# header keyword list, a curve, and an image.  It uses the XPA Tcl interface
# to control fv.  XPA is not distributed with fv, but can be obtained
# from the SAO/HEAD R&D group at 
#          http://hea-www.harvard.edu/RD/xpa/
# Be sure to build the XPA Tcl library.
#
# USAGE:
#
#   Execute this file using a command such as 'tclsh XPAScript.tcl' after
#   starting fv.  Or, start a TCL shell and do 'source XPAScript.tcl'.
#
# Another version of this script (TclScript.fv) uses fv's own Tcl macro
# interface, but it must be run by fv itself.
#
########################################################################

#
# Try loading the XPA Tcl Library
#

if { [catch {load [file join $env(LHEASOFT)/lib libtclxpa[info sharedlibextension]] tclxpa}] } {
   error "Cannot load XPA library"
}

#
# Set a variable here which points to the fits files to be opened.
# This is just to make it easier to specify files later.
#

set FitsDir ftp://heasarc.gsfc.nasa.gov/software/ftools/release/other/pdw

#
# Open 2 sample files
#

xpaset {} "fv" "open $FitsDir/ngc1316r.fit $FitsDir/rate.fit" \
      {} {} {} names errs 1

#
# Select one of the files and open a header window of extension #1
#

xpaset {} "fv" "select rate.fit" \
      {} {} {} names errs 1
xpaset {} "fv" "display header 1 " \
      {} {} {} names errs 1

#
# Plot a curve of Time vs Rate in POW and alter the graph's appearance
#

xpaset {} "fv" "display curve 1 time rate" \
      {} {} {} names errs 1
xpaset {} "fv" "pow bounds 770 -30 1070 300" \
      {} {} {} names errs 1
xpaset {} "fv" "pow curve pDisp No lDisp Yes lColor Blue" \
      {} {} {} names errs 1

#
# Select the other file and plot an image, setting its colormap to histogram
#

xpaset {} "fv" "select ngc1316r.fit" \
      {} {} {} names errs 1
xpaset {} "fv" "display image 0" \
      {} {} {} names errs 1
xpaset {} "fv" "pow colormap scale histo" \
      {} {} {} names errs 1
