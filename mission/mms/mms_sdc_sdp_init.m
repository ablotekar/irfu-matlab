function ENVIR = mms_sdc_sdp_init(scNumberStr)
% MMS_SDC_SDP_INIT reads initial environment and constants for MMS FIELDS processing
% 	[ENVIR, MMS_CONST] = MMS_SDC_SDP_INIT(scNumberStr) returns environment
%   variables and constants useful for MMS processing. Input argument 
%   should be the sc number (as a string), i.e. '1' for mms1 and '2' for 
%   mms2 etc. It also configures logging to "$LOG_PATH_ROOT/ mmsX/ sdp/
%   date_IRFU.log".
%
%       The struct ENVIR will contain the following:
%	  .DATA_PATH_ROOT       - Root dir of data files
%	  .CDF_BASE             - Root dir of CDF tools
%	  .DROPBOX_ROOT         - Root dir of our output files (temporary location)
%	  .LOG_PATH_ROOT        - Root dir of log files
% 	  .CAL_PATH_ROOT        - Root dir of calibration files
%
%	The struct MMS_CONST will contain the following:
%	  .Version.X		- Major Software version used.
%		  .Y		    - Major Calibration version used.
%		  .Z		    - File version (should perhaps be removed).
%	  .Bitmask.OnlyDCE  - Only DCE was found at these points in time. 
%
%	Example:
%		[ENVIR, MMS_CONST] = MMS_SDC_SDP_INIT('1');
%

narginchk(1,1); % SC number to ensure log is put in right place.

ENVIR = [];

ENVIR.CDF_BASE = getenv('CDF_BASE'); % Get path to CDF tools.
ENVIR.DATA_PATH_ROOT = getenv('DATA_PATH_ROOT'); % The final path of data.
ENVIR.LOG_PATH_ROOT = getenv('LOG_PATH_ROOT'); % Get path to logs.
ENVIR.DROPBOX_ROOT = getenv('DROPBOX_ROOT'); % Get path to output location,
% DROPBOX_ROOT is only the temporary location, file are then to be moved by
% other scripts to their final location when processing is done.
ENVIR.CAL_PATH_ROOT = getenv('CAL_PATH_ROOT'); % Get path to cal.

% Setup logging.
% Create a logfile at $LOG_PATH_ROOT / mmsX / sdp /
% named after current run day yyyymmdd and _IRFU.log. If this fails
% create it at $LOG_PATH_ROOT and include full date with seconds.
if( str2double(scNumberStr)>1 || str2double(scNumberStr)<4 )
    % Check to verify that output dir exists for logging. If not, create
    % it.
    if(~exist([ENVIR.LOG_PATH_ROOT, filesep, 'mms', scNumberStr, filesep, ...
            'sdp'], 'dir'))
        mkdir([ENVIR.LOG_PATH_ROOT, filesep, 'mms', scNumberStr], 'sdp');
    end
    irf.log('log_out', [ENVIR.LOG_PATH_ROOT, filesep, 'mms', ...
        scNumberStr, filesep, 'sdp', filesep, datestr(now,'yyyymmdd'),...
        '_IRFU.log']);
    mms_sdc_sdp_datamanager('init',str2double(scNumberStr))
    % Set log level
    irf.log('notice');
else
    irf.log('log_out', [ENVIR.LOG_PATH_ROOT, filesep, ...
        datestr(now,'yyyymmddTHHMMSS'), '_IRFU.log']);
    irf.log('debug');
    err_str = ['Matlab:MMS_SDC_SDP_INIT:InputArg scNumber incorrectly ',...
        'determined as: ',scNumberStr];
    irf.log('critical', err_str);
    error('Matlab:MMS_SDC_SDP_INIT',['MMS_SDC_SDP_INIT recieved an ', ...
        'unexpected sc number string. Input to processing should be ', ...
        'fullpath/filename.cdf according to MMS standard, mmsX_whatever',...
        ' where X = 1, 2, 3 or 4.']);
end

% Store information in the log file about which version of Matlab is used
% and which version of IRFU-MATLAB.
irf.log('notice', ['Matlab version used is: ' version]);
irf.log('notice', ['irfu-matlab version used is:', irf('version')]);
