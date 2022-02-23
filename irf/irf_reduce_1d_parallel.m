function PD =  irf_reduce_1d_parallel(obj,dim,x,varargin)
%% irf_reduce_parallel is same as PDIST.REDUCE Reduces (integrates) 3D distribution to 1D (line) only rearrage for prafor loop.
%   Example (1D):
%     f1D = reduce1D(iPDist, ''1D'',nDMPA,''vg'',vg1D,''nMC'',nMC);
%     irf_spectrogram(irf_panel('f1D'),f1D.specrec('velocity_1D'));
%
%
%   See more example uses in Example_MMS_reduced_ion_dist,
%   Example_MMS_reduced_ele_dist, and Example_MMS_reduced_ele_dist_2D
%
%   Options:
%     'nMC'    - number of Monte Carlo iterations used for integration,
%                for default number see IRF_INT_SPH_DIST
%     'base'   - set the base for the projection to cartesian 'cart'
%                (default) or polar 'pol' (only valid for 2D planes)
%     'vg'     - array with center values for the projection velocity
%                grid in [km/s], determined by instrument if omitted
%     'vg_edges' - array with edge values for the projection velocity
%                grid in [km/s]
%     'phig'   - array with center values for the projection
%                azimuthal angle in [rad]
%     'vint'   - set limits on the out-of-plane velocity to get
%                cut-like distribution in 2D or a cylindrical shell
%                in 1D in [km/s]
%     'aint'   - angular limit in out-of-plane direction to make
%                projection cut-like in 2D (only valid for 2D planes)
%     'scpot'  - sets all values below scpot to zero and changes the
%                energy correspondingly (only valid for electrons)
%     'lowerelim' - sets all values below lowerelim to zero, does not
%                change the energy. Can be single value, vector or
%                Tseries, for example 2*scpot
%     'weight' - how the number of MC iterations per bin is weighted,
%                can be 'none' (default), 'lin' or 'log'
%
%
%   The output is a PDist object with the reduced distribution where
%   'data' is the integrated phase space density and 'depend'
%   contains one (line) or two (plane) vectors of the velocity
%   centers. The units of the velocity is [km/s].
%
% The integration itself is performed in irf_int_sph_dist.m
%
% See also: IRF_INT_SPH_DIST, PDIST.PLOT_PLANE, PDIST.SPECREC,
% IRF_SPECTROGRAM
%
% To Do: (1) Merge with orignal class
%        (2) change documentation
%
%
% Author :    
% Optimize for parallel processing : Ajay Lotekar
%%
[~,args,nargs] = axescheck(varargin{:});
irf.log('warning','Please verify that you think the projection is done properly!');
if isempty(obj);
    irf.log('warning','Empty input.');
    return;
else,
    dist = obj;
end    % check the input return error if it is empty

% Check to what dimension the distribution is to be reduced
if any(strcmp(dim,{'1D','2D'}))
    dim = str2double(dim(1)); % input dim can either be '1D' or '2D'
else
    error('First input must be a string deciding projection type, either ''1D'' or ''2D''.')
end


if dim == 1 % 1D: projection to line
    
    % to chekc direction matrix and reshape as per the dist
    if isa(x,'TSeries')
        xphat_mat = x.resample(obj).norm.data;
    elseif isnumeric(x) && numel(size(x) == 3)
        xphat_mat = repmat(x,dist.length,1);
    elseif isnumeric(x) && all(numel(size(x) == [dist.length 3]))
        xphat_mat = x;
    end
    
    % make sure x is unit vector
    xphat_amplitude = sqrt(sum(xphat_mat.^2,2));
    if abs(mean(xphat_amplitude)-1) < 1e-2 && std(xphat_amplitude) > 1e-2 % make sure x are unit vectors,
        xphat_mat = xphat_mat./repmat(xphat_amplitude,1,3);
        irf.log('warning','|<x/|x|>-1| > 1e-2 or std(x/|x|) > 1e-2: x is recalculated as x = x/|x|.');
    end
elseif dim == 2 % 2D: projection to plane
    if isa(x,'TSeries') && isa(varargin{1},'TSeries')
        y = varargin{1}; varargin = varargin(2:end); % assume other coordinate for perpendicular plane is given after and in same format
        xphat_mat = x.resample(obj).norm.data;
        yphat_mat = y.resample(obj).norm.data;
    elseif isnumeric(x) && numel(size(x) == 3)
        y = varargin{1}; varargin = varargin(2:end); % assume other coordinate for perpendicular plane is given after and in same format
        xphat_mat = repmat(x,dist.length,1);
        yphat_mat = repmat(y,dist.length,1);
    elseif isnumeric(x) && all(numel(size(x) == [dist.length 3]))
        y = varargin{1}; varargin = varargin(2:end); % assume other coordinate for perpendicular plane is given after and in same format
        xphat_mat = x;
        yphat_mat = y;
    else
        error('Can''t recognize second vector for the projection plane, ''y'': PDistP.reduce(''2D'',x,y,...)')
    end
    
    % it's x and z that are used as input to irf_int_sph_dist
    % x and y are given, but might not be orthogonal
    % first make x and y unit vectors
    xphat_amplitude = sqrt(sum(xphat_mat.^2,2));
    yphat_amplitude = sqrt(sum(yphat_mat.^2,2));
    % These ifs are not really necessary, but could be there if one
    % wants to add some output saying that they were not put in
    % (inputted) as unit vectors. The definition of unit vectors is not
    % quite clear, due to tiny roundoff(?) errors
    if abs(mean(xphat_amplitude)-1) < 1e-2 && std(xphat_amplitude) > 1e-2 % make sure x are unit vectors,
        xphat_mat = xphat_mat./repmat(xphat_amplitude,1,3);
        irf.log('warning','|<x/|x|>-1| > 1e-2 or std(x/|x|) > 1e-2: x is recalculated as x = x/|x|.');
    end
    if abs(mean(yphat_amplitude)-1) < 1e-2 && std(yphat_amplitude) > 1e-2 % make sure y are unit vectors
        yphat_mat = yphat_mat./repmat(yphat_amplitude,1,3);
        irf.log('warning','|<y/|y|>-1| > 1e-2 or std(y/|y|) > 1e-2: y is recalculated as y = y/|y|.');
    end
    % make z orthogonal to x and y
    zphat_mat = cross(xphat_mat,yphat_mat,2);
    zphat_amplitude = sqrt(sum(zphat_mat.^2,2));
    zphat_mat = zphat_mat./repmat(zphat_amplitude,1,3);
    % make y orthogonal to z and x
    yphat_mat = cross(zphat_mat,xphat_mat,2);
    % check amplitude again, incase x and y were not orthogonal
    yphat_amplitude = sqrt(sum(yphat_mat.^2,2));
    if abs(mean(yphat_amplitude)-1) < 1e-2 && std(yphat_amplitude) > 1e-2  % make sure y are unit vectors
        yphat_mat = yphat_mat./repmat(yphat_amplitude,1,3);
        irf.log('warning','x and y were not orthogonal, y is recalculated as y = cross(cross(x,y),x)');
    end
    
    
    nargs = nargs - 1;
    args = args(2:end);
    
    % Set default projection grid, can be overriden by given input 'phig'
    nAzg = 32;
    dPhig = 2*pi/nAzg;
    phig = linspace(0,2*pi-dPhig,nAzg)+dPhig/2; % centers
end

% need to chek the 2D part in above loop

% make input distribution to SI units, s^3/m^6
dist = dist.convertto('s^3/m^6');

%% Check for input flags
% Default options and values
doTint = 0;
doLowerElim = 0;
nMC = 100; % number of Monte Carlo iterations
vint = [-Inf,Inf];
aint = [-180,180]; % azimuthal intherval
vgInput = 0;
vgInputEdges = 0;
weight = 'none';
correct4scpot = 0;
base = 'cart'; % coordinate base, cart or pol


if strcmp(dist.species,'electrons')
    isDes = 1
else,
    isDes = 0;
end

ancillary_data = {};
have_options = nargs > 1;


while have_options
    switch(lower(args{1}))
        case {'t','tint','time'} % time (undocumented, can be removed?)
            l = 2;
            tint = args{2};
            doTint = 1;
        case 'nmc' % number of Monte Carlo iterations
            l = 2;
            nMC = args{2};
            ancillary_data{end+1} = 'nMC';
            ancillary_data{end+1} = nMC;
        case 'vint' % limit on transverse velocity (like a cylinder) [km/s]
            l = 2;
            vint = args{2};
        case 'aint'
            l = 2;
            aint = args{2};
        case 'phig'
            l = 2;
            phig = args{2};
        case 'vg' % define velocity grid
            l = 2;
            vgInput = 1;
            vg = args{2}*1e3;
        case 'vg_edges'
            l = 2;
            vgInputEdges = 1;
            vg_edges = args{2}*1e3; % m/s
        case 'weight' % how data is weighted
            l = 2;
            weight = args{2};
            ancillary_data{end+1} = 'weight';
            ancillary_data{end+1} = weight;
        case 'scpot'
            l = 2;
            scpot = args{2};
            ancillary_data{end+1} = 'scpot';
            ancillary_data{end+1} = scpot;
            correct4scpot = 1;
        case 'lowerelim'
            l = 2;
            lowerelim = args{2};
            ancillary_data{end+1} = 'lowerelim';
            ancillary_data{end+1} = lowerelim;
            doLowerElim = 1;
            if isnumeric(lowerelim) && numel(lowerelim) == 1
                lowerelim = repmat(lowerelim,dist.length,1);
            elseif isnumeric(lowerelim) && numel(lowerelim) == dist.length
                lowerlim = lowerelim;
            elseif isa(lowerelim,'TSeries')
                lowerelim = lowerelim.resample(dist).data;
            else
                error(sprintf('Can not recognize input for flag ''%s'' ',args{1}))
            end
        case 'base' %
            l = 2;
            base = args{2};
    end
    args = args((l+1):end);
    if isempty(args), break, end
end

% while loop understood so far so good :)

% set vint ancillary data
ancillary_data{end+1} = 'vint';
ancillary_data{end+1} = vint;
ancillary_data{end+1} = 'vint_units';
ancillary_data{end+1} = 'km/s';

%% Get angles and velocities for spherical instrument grid, set projection
%  grid and perform projection
units = irf_units;
emat = double(dist.depend{1});

%---------------------------ATTENTION--------------------------------------
if doLowerElim
    lowerelim_mat = repmat(lowerelim, size(emat(1,:)));
else
    lowerelim_mat = [];  %# par for correction need to look into it
end
if correct4scpot
    scpot = scpot.tlim(dist.time).resample(dist.time);
    scpot_mat = repmat(scpot.data, size(emat(1,:)));
else
    scpot=[];      %# introduce due to par for loop
    scpot_mat =[]; %# introduced due to par for loop
end
%---------------------------END--------------------------------------------


if isDes == 1; M = units.me; else; M = units.mp; end
if doTint % get time indicies
    if length(tint) == 1 % single time
        it = interp1(dist.time.epochUnix,1:length(dist.time),tint.epochUnix,'nearest');
    else % time interval
        it1 = interp1(dist.time.epochUnix,1:length(dist.time),tint(1).epochUnix,'nearest');
        it2 = interp1(dist.time.epochUnix,1:length(dist.time),tint(2).epochUnix,'nearest');
        it = it1:it2;
    end
else % use entire PDistP
    it = 1:dist.length;
end
nt = length(it);
if ~nt % nt = 0
    error('Empty time array. Please verify the time(s) given.')
end


% try to make initialization and scPot correction outside time-loop

%---------------------------ATTENTION--------------------------------------
if not(any([vgInput,vgInputEdges])) % prepare a single grid outside the time-loop
    emax = dist.ancillary.energy(1,end)+dist.ancillary.delta_energy_plus(1,end);
    vmax = units.c*sqrt(1-(emax*units.e/(M*units.c^2)+1).^(-2));
    nv = 100;
    vgcart_noinput = linspace(-vmax,vmax,nv);
    irf.log('warning',sprintf('No velocity grid specified, using a default vg = linspace(-vmax,vmax,%g), with vmax = %g km/s.',nv,vmax*1e-3));
else
    vgcart_noinput=[]; %# Introduced due to par for loop
end
%---------------------------END--------------------------------------------

% loop to get projection
disp('Integrating distribution')
fprintf('it = %4.0f/%4.0f\n',0,nt) % display progress

%---------------------------ATTENTION--------------------------------------

if vgInputEdges
else
    vg_edges = [];
end


if doLowerElim
else
end
vg0 = vg;  %# Introduced due to par for loop
clear vg

% need to look into fg

Fg = zeros(length(it),length(vg0));  %# Introduced due to par for loop
vel = zeros(length(it),1);           %# Introduced due to par for loop
dens = zeros(length(it),1);          %# Introduced due to par for loop

%-------- taken from for loop
% elevation angle
th = double(dist.depend{3}); % polar angle in degrees
th = th-90; % elevation angle in degrees
th = th*pi/180; % in radi ans
%--------------------

%----- initializing the matrix ----
% if ~isempty(vg0)
%     all_vg = zeros(nt, length(vg0));
%     all_vg_edges = zeros(nt, length(vg0));
%
% end

%     all_vg =[];
%     all_vg_edges =[];


%-------------------
%---------------------------END--------------------------------------------


parfor i = 1:nt
    
    
    dispprog(i, nt) % replace with orignal loop ## if remove then can reduce 1 sec
    
    xphat = xphat_mat(i,:); % only use 1d arguments
    
    
    % 3d data matrix for time index it
    F3d = double(squeeze(double(dist.data(it(i),:,:,:)))); % s^3/m^6
    energy = emat(it(i),:);
    
    
    %perf_int(dist, xphat, F3d, energy, it, lowerelim_mat, scpot_mat
    
    if doLowerElim
        remove_extra_ind = 0; % for margin, remove extra energy channels
        ie_below_elim = find(abs(emat(it(i),:)-lowerelim_mat(it(i),:)) == min(abs(emat(it(i),:)-lowerelim_mat(it(i),:)))); % closest energy channel
        F3d(1:(max(ie_below_elim) + remove_extra_ind),:,:) = 0;
    end
    
    
    if correct4scpot
        if isfield(dist.ancillary,'delta_energy_minus') % remove all that satisfies E-Eminus<Vsc
            ie_below_scpot = find(emat(it(i),:)-dist.ancillary.delta_energy_minus(it(i),:)-scpot_mat(it(i),1)<0,1,'last');
            if 0 % disp energy channel that is removed, interferes with it = ... display
                disp(sprintf('Spacecraft potential = %g, Energy channel removed [E-Eminus,E,E+Eplus] = [%g,%g,%g]',...
                    scpot_mat(it(i),1),...
                    emat(it(i),ie_below_scpot)-dist.ancillary.delta_energy_minus(it(i),ie_below_scpot),...
                    emat(it(i),ie_below_scpot),...
                    emat(it(i),ie_below_scpot)+dist.ancillary.delta_energy_plus(it(i),ie_below_scpot)))
            end
        else
            ie_below_scpot = find(abs(emat(it(i),:)-scpot_mat(it(i),:)) == min(abs(emat(it(i),:)-scpot_mat(it(i),:)))); % closest energy channel
        end
        remove_extra_ind = 0; % for margin, remove extra energy channels
        F3d(1:(max(ie_below_scpot) + remove_extra_ind),:,:) = 0;
        %disp(sprintf('%8.1g ',energy))
        energy = energy-scpot_mat(it(i),:);
        %disp(sprintf('%8.1g ',energy))
        energy(energy<0) = 0;
        %disp(sprintf('%8.1g ',energy))
    end
    
    v = units.c*sqrt(1-(energy*units.e/(M*units.c^2)+1).^(-2)); % m/s
    
    % azimuthal angle
    if size(dist.depend{2},1)>1
        phi = double(dist.depend{2}(it(i),:)); % in degrees
    else % fast mode
        phi = double(dist.depend{2}); % in degrees
    end
    
    
    %phi = phi+180;
    %phi(phi>360) = phi(phi>360)-360;
    phi = phi-180;
    phi = phi*pi/180; % in radians
    
    
    % Set projection grid after the first distribution function
    % bin centers
    
    % need to add argument if vg is not initiate
    vg = ret_vg(vgInputEdges, vg_edges, vgInput, base, vgcart_noinput, v, vg0);
    
    
    % initiate projected f
    
    
    % perform projection  % removed 2d arguments
    
    % v, phi, th corresponds to the bins of F3d
    
    [all_v, all_vg_edg, tmpstF, tmpst_dens, tmpst_vel] = cal_all_vg(vgInputEdges, F3d,v,phi,th,vg,xphat,nMC,vint,aint,weight,vg_edges);
    all_vg(i,:) = all_v; % normally vg, but if vg_edges is used, vg is overriden
    
    all_vg_edges(i,:) = all_vg_edg;
    
    % fix for special cases % removed 2d arguments
    % dimension of projection, 1D if projection onto line, 2D if projection onto plane
    if dim == 1 || strcmpi(base,'cart')
        Fg(i,:,:) = tmpstF;
    end
    % set moments from reduced distribution (for debug)
    dens(i) = tmpst_dens;
    vel(i,:) = tmpst_vel;
    
    
end
%%
% Construct PDistP objects with reduced distribution
% vg is m/s, transform to km/s

PD = PDist(dist.time(it),Fg,'line (reduced)',all_vg*1e-3);
PD.ancillary.v_edges = all_vg_edges(1,:);
PD.species = dist.species;
PD.userData = dist.userData;
PD.ancillary.v_units = 'km/s';

% set units and projection directions
PD.units = 's/m^4';
PD.ancillary.projection_direction = xphat_mat(it,:);

while ~isempty(ancillary_data)
    PD.ancillary.(ancillary_data{1}) = ancillary_data{2};
    ancillary_data(1:2) = [];
end

if doLowerElim
    PD.ancillary.lowerelim = lowerelim_mat;
end
end

function dispprog(i, nt)
if mod(i,1) == 0,
    fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],i,nt);
end % display progress
end



function vg = ret_vg(vgInputEdges, vg_edges, vgInput, base, vgcart_noinput, v, vg)

if vgInputEdges % redefine vg (which is vg_center)
    vg = vg_edges(1:end-1) + 0.5*diff(vg_edges);
elseif vgInput
    vg = vg;
else % define from instrument velocity bins
    if strcmp(base,'cart')
        vg = vgcart_noinput; % maybe just bypass this and go directly through input vg_edges?
    else
        
        vg = [-fliplr(v),v]; % removed 2d argument if loop
    end
end
end


function [all_v, all_vg_edg, tmpstF, tmpst_dens, tmpst_vel]  = cal_all_vg(vgInputEdges, F3d,v,phi,th,vg,xphat,nMC,vint,aint,weight,vg_edges)

if vgInputEdges
    tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight,'vg_edges',vg_edges);
else
    tmpst = irf_int_sph_dist(F3d,v,phi,th,vg,'x',xphat,'nMC',nMC,'vzint',vint*1e3,'aint',aint,'weight',weight);
end
all_v = tmpst.v; % normally vg, but if vg_edges is used, vg is overriden
all_vg_edg(1,:) = tmpst.v_edges;
tmpstF = tmpst.F;

tmpst_dens = tmpst.dens;
tmpst_vel = tmpst.vel;
end