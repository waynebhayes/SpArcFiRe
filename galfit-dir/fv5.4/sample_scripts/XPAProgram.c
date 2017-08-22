/***********************************************************************
*
* Sample XPA/C program for fv 2.6
*
* This program will open 2 sample files (available online and distributed
* with the fv executables) and displays their contents in several forms...
* header keyword list, a curve, and an image.  It uses the XPA Tcl interface
* to control fv.  XPA is not distributed with fv, but can be obtained
* from the SAO/HEAD R&D group at 
*          http://hea-www.harvard.edu/RD/xpa/
*
* USAGE:
*
*   Compile this program using something like:
*         gcc -o xpa_c XPAProgram.c -I../../xpa-2.0 -L$FV/lib -lxpa
*   Then start fv and execute this program:
*         ./xpa_c
*
***********************************************************************/

#include <stdio.h>
#include <xpa.h>

#define NXPA 1

/*
 *   This function just cleans up any messages returned from fv
 *   after each XPA call and flags any errors.
 */
void PostProcessCall( int got, char **names, char **messages );

main()
{
   int  i, got, len, stat=0;
   char *buf,*fitsDir,cmd[255];
   char *names[NXPA];
   char *messages[NXPA];
   
   /*
    *  Set a variable here which points to the fits files to be opened.
    *  This is just to make it easier to specify files later.
    */
   
   fitsDir = "ftp://heasarc.gsfc.nasa.gov/software/ftools/release/other/pdw";

   /*
    *  Open 2 sample files
    */

   sprintf(cmd, "open %s/ngc1316r.fit %s/rate.fit", fitsDir, fitsDir);
   got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
   stat = PostProcessCall(got, names, messages);

   /*
    *  Select one of the files and open a header window of extension #1
    */

   if( !stat ) {
      sprintf(cmd, "select rate.fit");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

   if( !stat ) {
      sprintf(cmd, "display header 1");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

   /*
    *  Plot a curve of Time vs Rate in POW and alter the graph's appearance
    */

   if( !stat ) {
      sprintf(cmd, "display curve 1 time rate");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

   if( !stat ) {
      sprintf(cmd, "pow bounds 770 -30 1070 300");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

   if( !stat ) {
      sprintf(cmd, "pow curve pDisp No lDisp Yes lColor Blue");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

   /*
    *  Select the other file and plot an image, setting its colormap
    *  to histogram
    */

   if( !stat ) {
      sprintf(cmd, "select ngc1316r.fit");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

   if( !stat ) {
      sprintf(cmd, "display image 0");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

   if( !stat ) {
      sprintf(cmd, "pow colormap scale histo");
      got = XPASet(NULL, "fv", cmd, "", NULL, 0, names, messages, NXPA);
      stat = PostProcessCall(got, names, messages);
   }

}

int PostProcessCall( int got, char **names, char **messages )
{
   int i, status=0;

   for(i=0; i<got; i++){
      if( messages[i] != NULL ) {
         /* error processing */
         fprintf(stderr, "ERROR: %s (%s)\n", messages[i], names[i]);
         status = 1;
      }
      if( names[i]    ) free(names[i]);
      if( messages[i] ) free(messages[i]);
   }
   return status;
}
