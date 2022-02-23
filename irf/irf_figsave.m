function irf_figsave(varargin) %fname,resol,frm, ext)
%% Save figure in .jpg, .png, .eps and .pdf format 
%
% Input: fname = figure name
%        res = Resolution (default 300 dpi)
%        frem = frame size (default [12 6]
%        ext  = file extention (.png, .jpg, .pdf, .eps) (default .png)
%
% Example:
%       irf_figsave('figure01', 'res', 300, 'frem', [12 6], 'ext', '.png')
%       irf_figsave('figure01')
%%
args = varargin;
nargs = nargin;
orignal_args=args;

if nargs == 0 % show only help
    help irf_subplot;
    return
end

if nargs<2 || ~(mod(nargs,2)==0)
    error('IRFU_MATLAB:irf_figsave:InvalidNumberOfInputs','Incorrect number of input arguments')
end

%%%-------- Default Values ----------------
frm = [12 6];
resol = 300;
ext = '.png';

fname = args{1};     % file name

while ~isempty(args)
    x = args{1}; args(1) = [];
    switch lower(x)
        case {'res'}
            resol = args{1}; args(1) = [];
        case {'frem'} && length(args{1})==2
            frm = args{1}; args(1) = [];
        case {'ext'}
            ext= args{1}; args(1) = [];
    end
end


set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[frm(1) frm(2)]);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperPosition',[0 0 frm(1) frm(2)]);
set(gcf,'renderer','painters');
print(gcf, ['-d', ext] ,['-r', num2str(resol)],  fname)
end