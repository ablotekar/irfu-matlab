function irf_figsave(fname,resol,frm, ext)
%% Save figure in .jpg, .png, .eps and .pdf format 
%
%
% Example:
%       irf_figsave('figure01', 300, [12 6], '.png')
%
% Author : Ajay Lotekar

set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[frm(1) frm(2)]);
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperPosition',[0 0 frm(1) frm(2)]);
set(gcf,'renderer','painters');
print(gcf, ['-d', ext] ,['-r', num2str(resol)],  fname)
end
