function [out] = av_tsfilt(inp,fmin,fmax,Fs,order)
%function [out] = av_tsfilt(inp,fmin,fmax,[Fs],[order])
%
% inp, out   - column vectors
%             if inp has more than 1 column, 
%                                assume that the first column is time
%                                calculate frequency from the first time step
%                                assume that all time steps are the same length
% fmin,fmax  - filter frequencies
%               if fmin = 0 do lowpass filter
%               if fmax = 0 do highpass filter
% Fs         - sampling frequency if given as [] then Fs is calculated from time series 
% order      - the order of filter (elliptical IIR type filter is used) 
%              choose uneven order, 3 or 5 is OK. 
% 
% It is not better to have higher order filter. With orders above 10 for IIR filters like
% cheby1, ellip the result becomes wrong. For high pass filtering with very low 
% passband frequency (<.05) one can use [B2,A2] = cheby1(4,.3,.1,'high'); 
%
% Ex: def=av_tsfilt(de,0,.1,25,3); % lowpass filter E at .1Hz

warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_filt')

global AV_DEBUG; if isempty(AV_DEBUG), debug=0; else debug=AV_DEBUG;end

if ((nargin < 4) | (isempty(Fs))), 
 Fs=1/(inp(2,1)-inp(1,1));
 if debug == 1,disp(['Using sampling frequency ',num2str(Fs),' Hz']);end
end % estimate sampling frequency
if nargin > 4
    if debug == 1,disp(sprintf('You have specified %d-th filter order (use uneven order)',order));end
    n=order; % use this order for filters
end
fmin=fmin/(Fs/2);
fmax=fmax/(Fs/2);if (fmax > 1);fmax=1;end
out=inp;
Rp=.5;Rs=60;fact=1.1; % fact defines the width between stopband and passband
if fmin==0
  if fmax == 1, return;end
	if nargin < 5 
        [n wn]=ellipord(fmax,fmax*fact,Rp,Rs);
    end
	if debug == 1,disp(sprintf('using %d-th order lowpass filter',n));  end
	[B,A] = ellip(n,Rp,Rs,fmax);
elseif fmax ==0
    if nargin < 5
        [n wn]=ellipord(fmin,fmin*1.1,Rp,Rs);
    end 
	[B,A] = ellip(n,Rp,Rs,fmin,'high');
	if debug == 1,disp(sprintf('using %d-th highpass order filter',n));end
else
	%[n wn]=ellipord(fmax,fmax*1.1,Rp,Rs);
	%sprintf('using %d-th order ellip lowpass filter',n)
	%[B1,A1] = ellip(n,Rp,Rs,fmax);
	if nargin < 5
    	[n wn]=ellipord(fmax,fmax*1.3,Rp,Rs);
	end
	if debug == 1, disp(sprintf('using %d-th order ellip lowpass filter',n));end
	[B1,A1] = ellip(n,Rp,Rs,fmax);
	if nargin < 5
    	[n wn]=ellipord(fmin,fmin*.75,Rp,Rs);
	end
	if debug == 1,disp(sprintf('using	%d-th order ellip high pass filter',n));end
	[B2,A2] = ellip(n,Rp,Rs,fmin,'high');
end

[n m] = size(inp);
ini=1; % from which column start filtering
if m>1 % assume that first column is time
    ini=2;
end
if ((fmin ~= 0) & (fmax ~= 0))
	for i=ini:m
	out(:,i) = filtfilt(B1,A1,inp(:,i)); 
	out(:,i) = filtfilt(B2,A2,out(:,i)); 
	end
else
	for i=ini:m
	out(:,i) = filtfilt(B,A,inp(:,i)); 
	end
end

