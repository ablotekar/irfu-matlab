function res = find_diff(data,data0)
% find_diff(data,data0)
% returns the index of changing points
% Example:
% a=[0 0 1 1 1 1 0 0 0];
% >> find_diff(a')
% ans =
%
%     3
%     7
%
% $Id$
%
warning('caa:cleanup',...
'Function %s is deprecated and will be removed on May 1, 2004.\nUse %s instead',...
mfilename,'irf_find_diff')

a = size(data);
n_data = a(1);
if n_data < a(2), warning('check matrix orientation'), end

if nargin < 2
	data0 = data(1,:);
end

dd = data - ones(n_data,1)*data0;
der = abs(dd(2:end,:) - dd(1:end-1,:));

if a(2) > 1
	ii = find(der == 1) + 1;
	while 1
		ii_tocorr = find(ii > n_data);
		if isempty(ii_tocorr), break, end
		ii((ii_tocorr)) = ii((ii_tocorr)) - n_data;
	end
	res = rmDoubleIndex(ii);
else
	res = find(der == 1) + 1;
end

if dd(1)~=0
	if isempty(res)
		res = 1;
	else
		res(2:end+1) = res;
		res(1) = 1;
	end
end
