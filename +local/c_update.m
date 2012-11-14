function c_update
% LOCAL.C_UPDATE update index information in CAA directory
% 
% See also:
%	LOCAL.C_READ

% $Id$
% $Revision$  $Date$


dirCaa='/data/caa/CAA';
cd(dirCaa);
tmp=dir(dirCaa);
iDir = [tmp(:).isdir]; % find directories
dataSetArray = {d(iDir).name}';
dataSetArray(ismember(dataSetArray,{'.','..'})) = []; % remove '.' and '..'

for iDataSet=1:numel(dataSetArray)
	%% list files in data set directory
	dataSet=dataSetArray{iDataSet};
	listFiles=dir(dataSet);
	iDir = [listFiles(:).isdir]; %# returns logical vector
	listFiles(iDir)=[];
	%% read in file time intervals
	tmp=vertcat(listFiles.name);
	listFilesNames=[tmp tmp(:,end)]; % add one more column at the end
	listFilesNames(:,end)='=';       % changfe end column to character = (used as separator)
	listFilesNames=listFilesNames';
	tt=textscan(listFilesNames(:),'%*11s%4f%2f%2f_%2f%2f%2f_%4f%2f%2f_%2f%2f%2f%*s','delimiter','=');
	%% create index
	index.filename=[repmat([dataSet filesep],size(f,1),1) f];
	index.tstart=irf_time([tt{1} tt{2} tt{3} tt{4} tt{5} tt{6}],'vector2epoch');
	index.tend=irf_time([tt{7} tt{8} tt{9} tt{10} tt{11} tt{12}],'vector2epoch');
	eval(['index_' dataSet '=index;']);
	save('caa',['index_' datSet],'-append');
end

