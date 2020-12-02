%
% PSP_LOAD loads PSP data into TSeries
%
% PSP_LOAD(dataDir,datatype,dateStart,dateStop)
%
% PSP_LOAD(cdfFile,datatype) return variables into MATLAB base
%
% OUT = PSP_LOAD(..) output is TSeries if single variables is returned and
% cell array if several variables have to be returned.
%
% dataDir:  the directory with the cdf files
% datatype: 'mag_rtn','fgm_rtn'      - FGM, RTN coordinates
%           'mag_sc','fgm_sc'        - FGM, SC coordinates
%           'wf_dvdc'                - DBF Digital Fields Board Differential Voltage Waveform
%           'wf_scm'                 - DBF Digital Fields Board Search Coil Magnetometer Waveform
%           'ac_spec_dv12'           - DFB AC-coupled Differential Voltage V12 antenas
%           'ac_spec_dv34'           - DFB AC-coupled Differential Voltage V34 antenas
%           'ac_spec_v5'             - DFB AC-coupled V5 Antenna Voltage, Monopole Mode
%           'ac_spec_scmv'           - DFB AC-coupled SCM, Low Frequency, High Gain, v-component, Sensor coordinates
%           'ac_spec_scmu'           - DFB AC-coupled SCM, Low Frequency, High Gain, u-component, Sensor coordinates
%           'rfs_lfr'                - Radio Frequency Spectrometer, RFS, Low Frequency Reciever, LFR
%           'rfs_hfr'                - Radio Frequency Spectrometer, RFS, High Frequency Reciever, HFR
%           'sweap','spc'    - sweap proton moments
%           'ephem'          - ephemeris files
%           variable names as in the cdf file
% dateStart & dateStop: date vectors or strings for start and stop day
%            e.g. [yyyy mm dd] or 'yyyy/mm/dd' or  'yyyy mm dd' (read with datenum)
%
% Example:
%   psp_load('./','mag',[2020 01 26],[2020 02 01]); % assumes cdf files to be in the current directory
%   psp_load('psp_fld_l2_mag_RTN_20...cdf','mag');
%   rtnB = psp_load('psp_fld_l2_mag_RTN_20...cdf','mag');

function output = psp_load(dataDir,datatype,date_start,date_stop)

if nargout==0
  outputToBase = true;
else
  outputToBase = false;
end

switch datatype
  case {'mag_rtn', 'fgm_rtn'}
    filename= 'psp_fld_l2_mag_RTN';
    varnames = {'psp_fld_l2_mag_RTN'};
    varnamesout = {'rtnB'};
    
    hourtag={'00';'06';'12';'18'};
    
  case {'mag_sc', 'fgm_sc'}
    filename= 'psp_fld_l2_mag_SC';
    varnames = {'psp_fld_l2_mag_SC'};
    varnamesout = {'scB'};
    
    hourtag={'00';'06';'12';'18'};
    
  case {'wf_dvdc'}
    filename= 'psp_fld_l2_dfb_wf_dvdc';
    varnames = {...
      'psp_fld_l2_dfb_wf_dVdc_sensor';...
      'psp_fld_l2_dfb_wf_dVdc_sc'};
    varnamesout = {'wf_dvdc_sensor';'wf_dvdc_sc'};
    
    hourtag={'00';'06';'12';'18'};  
    
  case {'wf_scm'}
    filename= 'psp_fld_l2_dfb_wf_scm';
    varnames = {...
      'psp_fld_l2_dfb_wf_scm_hg_sensor';...
      'psp_fld_l2_dfb_wf_scm_hg_sc'};
    varnamesout = {...
      'wf_scm_sensor';'wf_scm_sc'};
    
    hourtag={'00';'06';'12';'18'};    
    
  case {'ac_spec_dv12'}
    filename = 'psp_fld_l2_dfb_ac_spec_dV12hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_dV12hg_frequency_bins';...          
      'psp_fld_l2_dfb_ac_spec_dV12hg'};
    varnamesout = {'ac_spec_dv12_freq_bins';'ac_spec_dv12_pw'};
    
    hourtag={''};  
    
  case {'ac_spec_dv34'}
    filename = 'psp_fld_l2_dfb_ac_spec_dV34hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_dV34hg_frequency_bins';...          
      'psp_fld_l2_dfb_ac_spec_dV34hg'};
    varnamesout = {'ac_spec_dv34_freq_bins';'ac_spec_dv34_pw'};
    
    hourtag={''};  
    
  case {'ac_spec_v5'}
    filename = 'psp_fld_l2_dfb_ac_spec_V5hg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_V5hg_frequency_bins';...          
      'psp_fld_l2_dfb_ac_spec_V5hg'};
    varnamesout = {'ac_spec_v5_freq_bins';'ac_spec_v5_pw'};
    
    hourtag={''};
    
  case {'ac_spec_scmv'}
    filename = 'psp_fld_l2_dfb_ac_spec_SCMvlfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_SCMvlfhg_frequency_bins';...          
      'psp_fld_l2_dfb_ac_spec_SCMvlfhg'};
    varnamesout = {'ac_spec_scmv_freq_bins';'ac_spec_scmv_pw'};
    
    hourtag={''};
    
  case {'ac_spec_scmu'}
    filename = 'psp_fld_l2_dfb_ac_spec_SCMulfhg';
    varnames = {...
      'psp_fld_l2_dfb_ac_spec_SCMulfhg_frequency_bins';...          
      'psp_fld_l2_dfb_ac_spec_SCMulfhg'};
    varnamesout = {'ac_spec_scmu_freq_bins';'ac_spec_scmu_pw'};
    
    hourtag={''}; 
    
    case {'rfs_lfr'}
    filename = 'psp_fld_l2_rfs_lfr';
    varnames = {...
      'psp_fld_l2_rfs_lfr_auto_averages_ch0_V1V2';...
      'frequency_lfr_auto_averages_ch0_V1V2';...
      'psp_fld_l2_rfs_lfr_auto_averages_ch1_V3V4';...
      'frequency_lfr_auto_averages_ch1_V3V4'};
    varnamesout = {'rfs_lfr_v1v2';'rfs_lfr_v1v2_freq';...
      'rfs_lfr_v3v4';'rfs_lfr_v3v4_freq'};
    
    hourtag={''};
  
  case {'rfs_hfr'}
    filename = 'psp_fld_l2_rfs_hfr';
    varnames = {...
      'psp_fld_l2_rfs_hfr_auto_averages_ch0_V1V2';...
      'frequency_hfr_auto_averages_ch0_V1V2';...
      'psp_fld_l2_rfs_hfr_auto_averages_ch1_V3V4';...
      'frequency_hfr_auto_averages_ch1_V3V4'};
    varnamesout = {'rfs_hfr_v1v2';'rfs_hfr_v1v2_freq';...
      'rfs_hfr_v3v4';'rfs_hfr_v3v4_freq'};
    
    hourtag={''};  
    
  case {'sweap', 'spc'}    
    filename = 'psp_swp_spc_l3i';
    
    varnames = {...
      'DQF';...          % data quality flag
      'general_flag';...
      'np_fit';'wp_fit';'vp_fit_SC';'vp_fit_RTN';...
      'np1_fit';'wp1_fit';'vp1_fit_SC';'vp1_fit_RTN';...
      'np_moment';'wp_moment';'vp_moment_SC';'vp_moment_RTN';...
      'sc_pos_HCI';'sc_vel_HCI'};
    varnamesout = {'DQF';'general_flag';...
      'np_fit';'wp_fit';'vp_fit_SC';'vp_fit_RTN';...
      'np1_fit';'wp1_fit';'vp1_fit_SC';'vp1_fit_RTN';...
      'np_moment';'wp_moment';'vp_moment_SC';'vp_moment_RTN';...
      'R_HCI';'V_HCI'};
    
    hourtag={''};
    
  case 'ephem'
    
    filename = 'spp_fld_l1_ephem_spp_rtn_';
    varnames = {'position';'velocity'};
    varnamesout = {'R_ephem';'V_ephem'};
    
    hourtag={''};
    
  otherwise
    if nargin==2
      % read in single variable from single file
      nFiles = 1;
      varName = datatype;
      filesToLoadTable = dataDir;
      [filename,varNameOut] = get_variable_filename(varName);
      if isempty(filename)
        error(['Filename not known for varname: ' varName '. Consider updating psp_load().']);
      elseif (filename ~= filesToLoadTable(1:length(filename)))
        error(['varname: ' varName ' not consistent with filename: ' filesToLoadTable]);
      end
      varnames = {varName};
      varnamesout = {varNameOut};
    else
      error('Data type not recognized!')
    end
end

if ~exist('filesToLoadTable','var')
  
  
  tStart = datenum(date_start);
  tStop  = datenum(date_stop);
  nDays  = tStop-tStart+1; % include first day
  
  datenumTable   = (tStart:tStop)';
  timestampTable = datestr(datenumTable,'yyyymmdd');
  
  % sanity check
  if size(timestampTable,1)~=nDays
    error('Number of days is odd.');
  end
  
  nHourtag = length(hourtag);
  nFiles = nDays*nHourtag;
  
  filesToLoadTable = char(zeros(nFiles,...
    size(timestampTable,2)+length(hourtag{1})));
  iFile = 1;
  for i = 1:nDays
    for i_hh=1:nHourtag
      filesToLoadTable(iFile,:) = [timestampTable(i,:) hourtag{i_hh}];
      iFile = iFile + 1;
    end
  end
  filesToLoadTable= strcat(dataDir,filesep,filename,'_',filesToLoadTable,'_v00','.cdf');
  
end

nVar      = length(varnames);
varData   = cell(nVar,1);
varMeta   = cell(nVar,1);
epochData = cell(nVar,1);

for iFile = 1:nFiles
  
  fileToLoad=strtrim(filesToLoadTable(iFile,:));
  irf.log('notice',['Loading: ' fileToLoad]);
  
  
  % version check
  fileFound = 0;
  nVersion = 20; % this is the highest version to check. Please update if ever needed.
  
  fileToLoad_vcheck = fileToLoad;
  for iVersion = 0:nVersion
      
      fileToLoad_vcheck(end-5:end-4) = num2str(iVersion,'%02d');
      
      if exist(fileToLoad_vcheck,'file')
          
          if fileFound == 1
                            irf.log('warning',['Version conflict, Replacing ''' fileToLoad ''' with ''' fileToLoad_vcheck ''''])
                            % give warning if two different versions of the same file exist, go tell someone, server should only keep the latest version.
          end
          
          fileToLoad = fileToLoad_vcheck;
          fileFound = 1;
      end
      
  end
  
      
  if fileFound
      
    pspobj=dataobj(fileToLoad);
    
    for iVar = 1:nVar
      
      varname    = varnames{iVar};
      varnameout = varnamesout{iVar};
      
      disp(['Loaded: ' fileToLoad ', Extracting: ' varname ', Output: ' varnameout ])
      
      var_data=get_ts(pspobj,varname);
      
      if isempty(var_data.data)
        error('Empty data!');
        
        % first data, allocate cellarray and get metadata
      elseif isempty(varData{iVar})
        varData{iVar}=cell(nFiles,1);
        epochData{iVar}=cell(nFiles,1);
        varMeta{iVar}=var_data;
      end
      
      varData{iVar}{iFile} = var_data.data;
      epochData{iVar}{iFile} =var_data.time.epoch;
      
    end
    
  end
  
  clear fileToLoad
  clear timestamp
  clear var_data
  clear pspobj
  
end

irf.log('notice','Done loading. Merging the data...');

if ~outputToBase
  output = cell(1,nVar);
end

% Fix the output 
for iOutputVar = 1:nVar
  
  varname_out = varnamesout{iOutputVar};
  indEmptyFiles = false(length(varData{iOutputVar}),1);
  for i=1:length(varData{iOutputVar})
    if isempty(varData{iOutputVar}{i})
      indEmptyFiles(i) = true;
    end
  end
  var_data = cell2mat(varData{iOutputVar}(~indEmptyFiles));
  epoch_data = cell2mat(epochData{iOutputVar}(~indEmptyFiles));
  varOut =  TSeries(EpochTT(epoch_data),var_data,...
    'metadata_from',varMeta{iOutputVar});
  
  if outputToBase
    assignin('base',varname_out,varOut);
  else
    output{iOutputVar} = varOut;
    if nVar==1
      output = output{1};
    end
  end
  
end

% GET_VARIABLE_FILENAME should return the filename for a given variable
% name and if variable has some shortened name. Currently only a few
% variable names are implemented
  function [fileName,shortVarOut] = get_variable_filename(varName)
    persistent pspVariables
    pspVariables=[...
      'psp_fld_l2_mag_RTN              psp_fld_l2_mag_RTN            rtnB       ';...
      'psp_fld_l2_mag_SC               psp_fld_l2_mag_SC             scB        ';...
      'psp_fld_l2_dfb_wf_V?dc          psp_fld_l2_dfb_wf_vdc                    ';... ? mark can be 1-5
      'psp_fld_l2_dfb_wf_scm_hg_sensor psp_fld_l2_dfb_wf_scm                    ';...
      'psp_fld_l2_dfb_wf_scm_hg_sc     psp_fld_l2_dfb_wf_scm                    ';...
      'psp_fld_l2_dfb_dc_spec_dV12hg   psp_fld_l2_dfb_dc_spec_dV12hg            ';...
      'psp_fld_l2_dfb_ac_spec_dV12hg   psp_fld_l2_dfb_ac_spec_dV12hg            ';...
      ];
    for iVariable = 1:size(pspVariables,1)
      varToTest = pspVariables(iVariable,1:length(varName));
      if all((varName == varToTest) | (varToTest == '?')) % allows to have any number in place of '?'
        fileName = strtrim(pspVariables(iVariable,33:58));
        shortVarOut = strtrim(pspVariables(iVariable,63:72));
        if isempty(shortVarOut)
          shortVarOut = varName;
        end
        irf.log('debug',['Variable: ' varName '  File Name: ' fileName]);
      end
    end
  end

end
