#!/bin/bash

# This is used for starting the MATLAB process of irfu-matlab code on SDC for the MMS mission.
# Author: T. Nilsson, IRFU
# Date: 2014/02/28
#
# Usage: place script in the same folder as has irfu-matlab as a subfolder, then run
#  "./script.sh <mmsX_dce_filename> <mmsX_dcv_filename> <mmsX_101_filename>", with the following
#  input arguments: (order is irrelevant as long as the OptionalDataDescriptor is "_dce_", "_dcv_", or "_101_").
#    <mmsX_***_dce_filename.cdf> = Filename of DC E data to be processed for SITL. Including path and extension.
#    <mmsX_***_dcv_filename.cdf> = Filename of DC V data to be processed for SITL. Including path and extention.
#    <mmsX_***_101_filename.cdf> = Filename of HK 101 data (with sunpulse) to be processed for SITL. Including path and extention.
#  output files created:
#    <mmsX_***_sitl_***_dce2d_yyyymmddHHMMSS_vX.Y.Z.cdf>         = File placed in $DROPBOX_ROOT
#    <DATE_IRFU.log>                                             = Logfile of run, placed in $LOG_PATH_ROOT/mmsX/sdp/.
#    <mmsX_***_sitl_***_dce2d_yyyymmddHHMMSS_vX.Y.Z_runTime.log> = File to identify output file, as per e-mail of 2013/11/22. Also placed in $LOG_PATH_ROOT/mmsX/sdp/.
#
#    return code 0 if ok.
#    return code 199 if built in Matlab code error.
#    return code 100-198 if MMS Matlab specific code error. <<-- FIXME: These are still to be defined and specified.
#
# Note: The script assumes it is located in the folder which has irfu-matlab as a subfolder.

# make sure that the correct number of arguments are provided
if [ ${#} -lt 2 ] || [ ${#} -gt 3 ] ; then
  exit 166  # SDC-defined error code for "incorrect usage"
fi

# For debug: display input arguments to terminal
#echo $1 # First argument
#echo $2 # Second argument
#echo $3 # Third argument

# SET ENVIRONMENT MATLABPATH by finding all dir and subdir of 'irfu-matlab' with full path.
#find `pwd` -type d \( -name '@*' -o -name '+*' -o -name '.git' \) -prune -o \( -path "*irfu-matlab*" -type d \) -printf %p:
export MATLABPATH="$(find `pwd` -type d \( -name '@*' -o -name '+*' -o -name '.git' \) -prune -o \( -path "*irfu-matlab*" -type d \) -printf %p:)$MATLABPATH"

if [ ${#} -eq 2 ] ;  then
# RUN THIS IF ONLY ONE FILE EXISTS (DCE)
   matlab -nodesktop -nosplash -nodisplay -nojvm -r "try, mms_sdc_sdp_proc('sitl','$1','$2'), catch err, if(strcmp(err.identifier,'MATLAB:SDCcode')) irf.log('critical',['Bash error catch worked: ', err.identifier, '. With message: ', err.message]); exit(str2num(err.message)); else irf.log('critical',['Bash error catch: ', err.identifier, '. With message: ', err.message]); exit(199); end; end, exit(0)"
elif [ ${#} -eq 3 ] ; then
# RUN THIS IF TWO FILES EXISTS (DCE and DCV)
   matlab -nodesktop -nosplash -nodisplay -nojvm -r "try, mms_sdc_sdp_proc('sitl','$1','$2','$3'), catch err, if(strcmp(err.identifier,'MATLAB:SDCcode')) irf.log('critical',['Bash error catch worked: ', err.identifier, '. With message: ', err.message]); exit(str2num(err.message)); else irf.log('critical',['Bash error catch: ', err.identifier, '. With message: ', err.message]); exit(199); end; end, exit(0)"
fi

# For debug: display message 'all done'.
#echo $? # Return code 0 if ok, 100-199 if error occured.
#echo 'Back in bash, all done.'
