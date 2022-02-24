function  PVI = irf_compute_pvi(B, avhr, Ntau, dint)
%%
% This function compute the PVI 
% 
% PVI = irf_compute_pvi(B, avhr, Ntau, dint)
% 
% Input: B = B field structure B.epoch and B.data
%       avhr = avarage period for sigma (in hr)
%       Ntau = Number of points over PVI wants compute
%       dint = data intervel 
% 
% Output : PVI.epoch = PVI time 
%         PVI.PVI = B increment along different component (PVI.PVI.x ...)
%         PVI.sigma = different compoent of sigma
%         PVI.normPVI = normalize PVI
%         PVI.indexPVI = PVI index
%
% Author : Ajay Lotekar
%%


dt=min(diff(B.epochUnix));   % Minimum time step

PVIx=B.data((Ntau+1):end, 1)-B.data(1:end-Ntau, 1);     % dBx
PVIy=B.data((Ntau+1):end, 1)-B.data(1:end-Ntau, 1);     % dBy
PVIz=B.data((Ntau+1):end, 1)-B.data(1:end-Ntau, 1);     % dBz
PVI=sqrt(PVIx.^2+PVIy.^2+PVIz.^2);            % Dinominator in above equation
time_PVI=B.epochUnix(1:end-Ntau);


tmin=min(B.epochUnix); tmax=max(B.epochUnix);    % min and max To values
%time step
Jo={};
avhr = 2;
dT=2*60*60; % avarage sampling time min*60
for i=1:floor(dint/avhr)
    [I,Ja]=find(time_PVI'>=tmin+(i-1)*dT); [I,Jb]=find(time_PVI'<=tmin+i*dT); Jo{i}=intersect(Ja,Jb);
end

norm_PVIx=zeros(size(PVIx));
norm_PVIy=zeros(size(PVIy));
norm_PVIz=zeros(size(PVIz));

sigma_x=[]; sigma_y=[]; sigma_z=[];
for i=1:floor(dint/avhr)
    range=Jo{i};
    time_sigma=mean(time_PVI(range));
    sigma_x=[sigma_x sqrt(mean(PVIx(range).^2))];
    sigma_y=[sigma_y sqrt(mean(PVIy(range).^2))];
    sigma_z=[sigma_z sqrt(mean(PVIz(range).^2))];
    
    norm_PVIx(range)=PVIx(range)./sigma_x(end);
    norm_PVIy(range)=PVIy(range)./sigma_y(end);
    norm_PVIz(range)=PVIz(range)./sigma_z(end);
end

norm_PVI=sqrt(norm_PVIx.^2+norm_PVIy.^2+norm_PVIz.^2);

PVI = [];
PVI.epoch = time_PVI;
PVI.PVI.x=PVIx;
PVI.PVI.y=PVIy;
PVI.PVI.z=PVIz;

PVI.sigma.x = sigma_x;
PVI.sigma.y = sigma_y;
PVI.sigma.z = sigma_z;

PVI.normPVI.x = norm_PVIx;
PVI.normPVI.y = norm_PVIy;
PVI.normPVI.z = norm_PVIz;

PVI.indexPVI = norm_PVI;

end