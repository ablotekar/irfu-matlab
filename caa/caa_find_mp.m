function [t_mp_out,t_mp_in] = caa_find_mp(start_time, dt, cl_id)
%CAA_FIND_MP  find model magnetopause crossings
%
% [t_mp_out,t_mp_in] = caa_find_mp(start_time, dt, cl_id)
%
% See also IRF_SHUE_MP
%
% $Id$

% Copyright 2006 Yuri Khotyaintsev

if dt>toepoch([1996 01 01 00 00 00])
	% et is given
	if dt< start_time, error('STOP_TIME must be larger then START_TIME)'), end
	dt = dt - start_time;
end

t_mp_out = []; t_mp_in = [];

R_E = 6378;
ACE_X_POS = 222*R_E;	% ACE X postition
ACE_VX_DEF = 480;		% Default solar wind speed
ACE_DT_DEF = ACE_X_POS/ACE_VX_DEF;
ACE_N_DEF = 1;			% Defaultsolar wind density
ACE_BZ_DEF = 0;			% Default IMF Bz

irf_log('proc',['orbit : ' epoch2iso(start_time,1) ' -- ' ...
		epoch2iso(start_time+dt,1)])

data = getData(ClusterDB, start_time, dt, cl_id, 'r', 'nosave');
if isempty(data), error('cannot fetch position'), end

R = data{2};
R = R(R(:,1)>0,:); % we probably cross the MP only for positive X
R = R(irf_abs(R,1)>7*R_E,:); % we probably cross the MP only R > 7 R_E

if isempty(R)
	irf_log('proc','tail season')
	return
end

start_time = R(1,1);
dt = R(end,1) -R(1,1);

irf_log('proc',['X>0, R>7R_E: ' epoch2iso(start_time,1) ' -- ' ...
		epoch2iso(start_time+dt,1)])

% Fetch ACE data
ISTP_PATH = '/data/istp';
ace_B = irf_istp_get(ISTP_PATH, start_time -60*60, dt +120*60, 'ace', 'b');
ace_V = irf_istp_get(ISTP_PATH, start_time -60*60, dt +120*60, 'ace', 'v');
ace_N = irf_istp_get(ISTP_PATH, start_time -60*60, dt +120*60, 'ace', 'n');

% Create new timeline with 30 min step
st_a = fromepoch(start_time);
st = toepoch([st_a(1:4) fix(st_a(5)/30)*30 00]);
dt = ceil((start_time +dt -st)/1800)*1800;

irf_log('proc',['subint: ' epoch2iso(start_time,1) ' -- ' ...
		epoch2iso(start_time+dt,1)])

%v_ttt = []; b_ttt = []; n_ttt = [];
r_prev = [];
for t=st:1800:st+dt
	%irf_log('proc',['time: ' epoch2iso(t,1)])
	
	% ACE time shift
	if isempty(ace_V), ace_dt = ACE_DT_DEF;
	else
		v_tmp = linear_solve(ace_V, t, ACE_DT_DEF);
		%irf_log('proc',['ace_v_tmp: ' num2str(round(v_tmp)) ' km/s'])
		dt_ace = ACE_X_POS/v_tmp;
	end
	%irf_log('proc',['ace_dt   : ' num2str(round(dt_ace/60)) ' min'])
	
	if isempty(ace_V), vx_tmp = ACE_VX_DEF;
	else vx_tmp = linear_solve(ace_V, t, dt_ace);
	end
	%irf_log('proc',['ace_vx_tmp: ' num2str(vx_tmp,'%.2f') ' km/s'])
	%v_ttt = [v_ttt; t-dt_ace vx_tmp];
	
	if isempty(ace_V), n_tmp = ACE_N_DEF;
	else n_tmp = linear_solve(ace_N, t, dt_ace);
	end
	%irf_log('proc',['ace_nn_tmp: ' num2str(n_tmp,'%.2f') ' cc'])
	%n_ttt = [n_ttt; t-dt_ace n_tmp];
	
	if isempty(ace_V), bz_tmp = ACE_BZ_DEF;
	else bz_tmp = linear_solve(ace_B(:,[1 4]), t, dt_ace);
	end
	%irf_log('proc',['ace_vx_tmp: ' num2str(bz_tmp,'%.2f') ' nT'])
	%b_ttt = [b_ttt; t-dt_ace bz_tmp];
	
	r_tmp = linear_solve(R, t, 0);
	r_gsm = irf_gse2gsm([t r_tmp]);
	r_gsm(2:4) = r_gsm(2:4)/R_E;
	r_mp = irf_shue_mp(r_gsm, bz_tmp, nv2press(n_tmp,vx_tmp^2)); 
	%irf_log('proc',['r: ' num2str(r_gsm(2:4),'%.2f %.2f %.2f') ...
	%		' mp:' num2str(r_mp,'%.2f') ' Re'])
			
	if isempty(t_mp_out) && ~isempty(r_prev) && (r_prev>0) && (r_mp<0)
		t_mp_out = t -1800;
		irf_log('proc',['FOUND OUTBOUND : ' epoch2iso(t_mp_out,1)])
	end
	if ~isempty(r_prev) && (r_prev<0) && (r_mp>0)
		t_mp_in = t;
		irf_log('proc',['FOUND INBOUND  : ' epoch2iso(t_mp_in,1)])
	end
	r_prev = r_mp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Help functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = linear_solve(f, t, delta_t)
% help function for linear interpolation

t = t -delta_t;
p1 = f( (f(:,1) > t -3600) & (f(:,1) <= t), : ); % data has 1 hour resolutuion
p2 = f( (f(:,1) < t +3600) & (f(:,1) > t), : );
if isempty(p1) && isempty(p2), error('bizzare no data'), end
if ~isempty(p1), p1 = p1(1,:); end
if ~isempty(p2)
	p2 = p2(end,:);
	if isempty(p1), y = p2(:,2:end);
	else
		a_tmp = (p1(:,2:end) -p2(:,2:end))/(p1(:,1) -p2(:,1));
		y = a_tmp*t + 0.5*(p1(:,2:end) +p2(:,2:end) ...
				-a_tmp*(p1(:,1) +p2(:,1)));
	end
else y = p1(:,2:end);
end
return

function res = nv2press(n,v2)
%function res = nv2press(n,v2)
%
% Calculate plasma dynamic pressure in nPa
% n in 1/cc
% v^2 in [km/s]^2

n = n(:);
v2 = v2(:);
% p=nmv^2 ;-)
res = 1.6726*1e-6*v2.*n;
return
