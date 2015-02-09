function phase = mms_sdc_sdp_comp_phase
%MMS_SDC_SDP_COMP_PHASE  compute phase for datamanager
%
%  phase = mms_sdc_sdp_comp_phase()

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end
phase = MMS_CONST.Error;

procId = mms_sdc_sdp_datamanager('procId');
switch procId
  case {MMS_CONST.SDCProc.scpot,MMS_CONST.SDCProc.sitl, MMS_CONST.SDCProc.ql}
    hk_101 = mms_sdc_sdp_datamanager('hk_101');
    if isnumeric(hk_101) && numel(hk_101)==1 && hk_101==MMS_CONST.Error,
      irf.log('warning','Bad hk_101 input'); return
    end
    dce = mms_sdc_sdp_datamanager('dce');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad dce input'); return
    end
    [dcephase, dcephase_flag] = mms_sdc_sdp_phase_2(hk_101, dce.time);
    phase = struct('data',dcephase,'bitmask',dcephase_flag);
  case MMS_CONST.SDCProc.l2pre
    %XXX: this needs to be replaced with definitive attitude!!!
    hk_101 = mms_sdc_sdp_datamanager('hk_101');
    if isnumeric(hk_101) && numel(hk_101)==1 && hk_101==MMS_CONST.Error,
      irf.log('warning','Bad hk_101 input'); return
    end
    dce = mms_sdc_sdp_datamanager('dce');
    if isnumeric(dce) && numel(dce)==1 && dce==MMS_CONST.Error,
      irf.log('warning','Bad dce input'); return
    end

%   % DEFATT
%   defatt = mms_sdc_sdp_datamanager('defatt');
%   if isnumeric(defatt) && numel(defatt)==1 && defatt==MMS_CONST.Error,
%     irf.log('warning','Bad DEFATT input'); return
%   end
%
%   % Every time differance is less than 180 degrees assume it wrap'ed once
%   % and add 360 degrees, cumulatively.
%   % Assumption, DEFATT datapoints at least every 10 second
%   % (nominal spinrate once every 20 seconds) and all s/c spinning in +Z
%   % direction.
%   phase_jumps = [0; diff(defatt.zphase)<-180]*360;
%   defatt.zphase = defatt.zphase + cumsum(phase_jumps);

% FIXME phase calculation primarly from polyfit() and polyval if fit was
% good, otherwise use interp1(). S gives normal residual (ie indicator of
% quality), mu is scaling and centering.
% [p, S, mu] = polyfit(double(defatt.time), defatt.zphase, 1);
% if (S.normr <= 5)
%  % "Good" fit, ie almost constant spin rate. Use polyval.
%  [dcephase, delta] = polyval(p, double(dce.time), S, mu);
%  % FIXME Set bitmask to indicate phase from DEFATT and polyval, possibly
%  % use delta for indicating quality of phase.
% else
%   % "Bad" fit, somewhat changing spin rate. Use table lookup via
%   % interp1(), linear interpolation with extrapolation.
%   dcephase = interp1(double(defatt.time), defatt.zphase, double(dce.time), 'linear', 'extrap');
%   % FIXME Set bitmask to indicate phase from DEFATT but with interp1.
% end
% % Then wrap to interval [0 - 360) degrees.
% dcephase = mod(phase, 360); % interval [0 to 360)

    [dcephase, dcephase_flag] = mms_sdc_sdp_phase_2(hk_101, dce.time);
    phase = struct('data',dcephase,'bitmask',dcephase_flag);


  case MMS_CONST.Error
    errStr = 'mms_sdc_sdp_datamanager not properly initialized';
    irf.log('critical',errStr), error(errStr)
  otherwise
    errStr = 'unrecognized procId';
    irf.log('critical', errStr); error(errStr)
end
