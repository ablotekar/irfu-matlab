function [DCE_SRF,PSP,ScPot,codeVerStr,matVerStr] = vdccal(VDC,EDC)
%SOLO.VDCCAL  Calibrate VDC to get DC E and PSP
%
%    [DCE_SRF,PSP,ScPot,codeVerStr,matVerStr] = solo.vdccal(VDC,EDC)
%
% Inputs: VDC,EDC from L2 CWF files
%
% Outputs:
%   DCE_SRF    - DC electric field in SRF (Ex=0)
%   PSP        - Probe-to-spacecraft potential
%   ScPot      - Spacecraft potential (PRELIMINARY PROXY)
%   codeVerStr - Date format version string for function itself. Used by BICAS.
%   matVerStr  - Date format version string for .mat file. Used by BICAS.
%                (Not yet used.)
%
% Loads .mat file produced by solo.correlate_probes_batch (script)
% NOTE: .mat needs to be updated before processing the new data.
%
% NOTE: This function is used by BICAS for producing official datasets.

a = load('d23K123_20210129');



%===========================================================================
% Date strings that represent the version of calibration. These strings are
% used by BICAS to set a CDF global attribute in official datasets for
% traceability.
% --
% Version of the function (not .mat file).
% NOTE: This value is meant to be be updated by hand, not by an automatic
% timestamp, so that a constant value represents the same function.
%===========================================================================
codeVerStr = '2021-04-08T14:00:00';
% Version of the .mat file. Not yet used. Meant to be read from .mat file.
matVerStr  = [];



Gamma0 = a.Gamma0;
Gamma1 = a.Gamma1;
cc = a.CC;



%=============================================================================
% Find data points/CDF records for which only V1_DC is available
%
% NOTE: solo.vdccal() only uses DC data (not AC). For mux=0, there are two
% cases:
% (1) all single probes (DC) available
% (2) only probe 1 (DC) available.
% NOTE: Only works for mux=0,2,3,4 (not mux=1).
% NOTE: Ignores EDC argument.
%==============================================================================
bSingleProbe = isnan(VDC.y.data) & isnan(VDC.z.data);



d23R  = a.d23.resample(VDC);
K123R = a.K123.resample(VDC);
Gamma1R = Gamma1.resample(VDC);
Gamma0R = Gamma0.resample(VDC);
ccR = cc.resample(VDC);

V2corr = double(VDC.y.data) -double(d23R.data); %Remove potential offset between 2,3
V23_corr = (V2corr+double(VDC.z.data))/2; %(V2corr+V3)/2
V2cmr = double(V2corr)-(Gamma0R.data+V23_corr.*Gamma1R.data)/2; %Remove common mode from V2.
V3cmr = double(VDC.z.data)+(Gamma0R.data+V23_corr.*Gamma1R.data)/2;

V23 = (V2cmr + V3cmr)/2; % (V2cmr + V3) /2

V23corr = (V23.*K123R.data(:,1) + K123R.data(:,2)); %Correcting V23 to V1


PSP = irf.ts_scalar(VDC.time,(V23corr + double(VDC.x.data))/2); %Compute PSP from corrected quantities.

% Use alternate, simpler "calculation" for single-probe data.
PSP.data(bSingleProbe) = VDC.x.data(bSingleProbe);

PSP.units = 'V';



PLASMA_POT = 1.5; SHORT_FACTOR = 2.5; % XXX: these are just ad hoc numbers.

ScPot = irf.ts_scalar(VDC.time, -(PSP.data-PLASMA_POT)*SHORT_FACTOR);
ScPot.units = PSP.units;

% Ey_SRF = V3 - V2, 6.99 - 1/2 of distance between the antennas
V_delta23_corr = V2cmr-V3cmr; %Fixed V2-V3.
Ey_SRF = -V_delta23_corr*1e3/6.99;

% Ez_SRF = V23corr - V1
% Here we use the antenna length of 11.2 m, which correponds to 
% the distance between the center of ANT1 and a symmetric antenna on the 
% other side having voltage V23 corr.
Ez_SRF = V23corr - double(VDC.x.data);	
Ez_SRF = Ez_SRF*1e3/11.2;

DCE_SRF = irf.ts_vec_xyz(VDC.time,[Ey_SRF*0 Ey_SRF Ez_SRF]);
DCE_SRF.units = 'mV/m';
DCE_SRF.coordinateSystem = 'SRF';
