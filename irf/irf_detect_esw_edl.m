function SORTM = irf_detect_esw_edl(varargin)
% DETECT_ESW_DEL this function detect bipolar and monopolar pulses from
% given signal
%
% SORTM = irf_detect_esw_edl(Em, tm, figpath, './figure/mms1', 'crESW', 0.8, 'crEDL', 0.8)
% Detect the bipolar and monopolar pulses in given signal.
%
%Input
% Em = Electric field (E parallel )
% tm = time
% figpath = path to save sorted structure.
% crESW = cross correlation between theory and observation bipolar pulse 
% crEDL = cross correlation between theory and observation monopolar pulse 
% Default value of
% crESW = 0.8;
% crEDL = 0.8;
%
%
% Output
% SORTM = Contain the inforation of the onset and left and righ bound and
%         width to take if want to replot sorted structure
% SORTM(:,1) = left bound (sec)
% SORTM(:,2) = right bound (sec)
% SORTM(:,3) = onset time (sec)
% SORTM(:,3) = width  (Number of grid points) to convert in sec multiply by
%              resolution of data

%
% Examples
% SORTM = irf_detect_esw_edl(Em, tm, 'figpath', './figure/mms1')
%
% To Do : chaneg it for structure
%
%
% Author: Ajay Lotekar
%%
SORTM = []; % Initializing output varible

args = varargin;
nargs = nargin;
orignal_args=args;

if nargs == 0 % show only help
    help irf_detect_esw_edl;
    return
end

if nargs<2 || ~(mod(nargs,2)==0)
    error('Incorrect number of input arguments')
end

% ---- default values --------------
figpath=[];
ampL = 1;    % Threshould amplitude of E field (mV/m)
crESW = 0.8; % CCR with ESW
crEDL = crESW; % CCR with EDL
% ----------------------------------
Epf = args{1}; args(1) = [];
tm = args{1}; args(1) = [];

while ~isempty(args)
    x = args{1}; args(1) = [];
    switch lower(x)
        case lower({'figpath'})
            figpath = args{1}; args(1) = [];
        case lower({'crESW'})
            crESW = args{1}; args(1) = [];
        case lower({'crEDL'})
            crEDL = args{1}; args(1) = [];
    end
end

if isempty(figpath)
    mkdir('./figure/sorted')
    figpath = './figure/sorted';
else
    mkdir(figpath)
end
%%

Epa = abs(Epf);  % absolute value for the detecting positve as well as negative peakes
Gm = (1:length(Epa)); % Grid point number
dt = min(diff(tm));

%%

multWD = 5;  % Multiply l0 by factor
Jold=9999;
Rbold = 1;
%------ Memory Initialization ---------------
dif(1:7)=0;     % part of peak detection algorithum
KJ = 1;         % some flag
%------------- for ploting --------------------------
fsiz = 12;
lw   = 2;

figure('Visible', 'off')

%  -------------------------------------------
for J=4:length(Epa)-4
    if Epa(J) >= ampL
        for IJ =-3:3
            dif(IJ+4) = Epa(J+IJ+1)-Epa(J+IJ);
        end
        %-----------------Find peak conditions ---------------------
        cond1 =  (dif(2)>0) && (dif(3)>0); % (dif(1)>0) &&
        cond2 = (dif(5)<0) && (dif(6)<0) ; % && (dif(7)<0)
        cond3 = cond1 && cond2;
        %------------------------------------------------------
        if cond3
            if abs(Jold-J) > 1
                
                Jold = J;
                for I=J:length(Epa)
                    if (Epa(I) <= 0.2*Epa(J))
                        idx1 = abs(I-J);
                        break
                    end
                end
                l0= idx1;
                if l0 <5
                    wd = multWD*10;
                else
                    wd = multWD*round(l0,-1);
                end
                
                if (l0 > 0.0) && (J-wd > 0) && (J+wd < length(Epf))
                     %disp(J+wd)
                    Es = Epf(J-wd:J+wd);
                    Gms = Gm(J-wd:J+wd);
                    ts = tm(Gms);
                    dx = -wd:1:wd;
                    
                    ph = exp(-((dx)/l0).^2);
                    Et = diff(ph)./diff(dx);
                    a0 = max(Et);
                    Et = Epa(J)*Et/a0;
                    
                    Es(length(Es))=[];
                    Gms(length(Es))=[];
                    dx(length(dx))=[];
                    %------------------------------------------------------
                    [mxa, a] = max(Et);
                    [mib, b] = min(Et);
                    mimx = [a, b];
                    l=min(mimx); % left bound
                    r=max(mimx); % right bound
                    Etct = Et(l:r);
                    [c, d] = min(abs(Etct));
                    t0c =l+d; % center
                    %------------------------------------------------------
                    for ic = wd:length(Et)
                        if ph(ic)<=0.01
                            ind1 = ic-wd;
                            break
                        end
                    end
                    lb = t0c-ind1;%min([ind1, ind2]);
                    rb = t0c+ind1;%max([ind1, ind2]);
                    
                    [corC,lagsT] = xcov(Es, Et,'coeff');
                    [m_cr, ncr] = max(abs(corC));
                    sng = corC(ncr)/abs(corC(ncr));
                    Et = sng*Et;
                    ampH = max(Et)+0.1*max(Et);
                    Etdl = exp(-((dx)/(0.5*l0)).^2);
                    Etdl = Epf(J)*Etdl;
                    [corC,lagsL] = xcov(Es, Etdl,'coeff');
                    [m_crdl, ncrdl] = max(abs(corC));
                    %------------------------------------------------------
                    if (lb+lagsT(ncr)) <= 0 || (rb+lagsT(ncr)) >=length(Es)+1
                        continue;
                    end
                    %------------------------------------------------------
                    Emx = max(Es(lb+lagsT(ncr):rb+lagsT(ncr)));
                    Emn = min(Es(lb+lagsT(ncr):rb+lagsT(ncr)));
                    ampM = [Emx, abs(Emn)];
                    ampR = min(ampM)/max(ampM);
                    Lb = Gms(lb+lagsT(ncr));
                    Rb = Gms(rb+lagsT(ncr));
                    T0 = Gms(t0c+lagsT(ncr));
                    %------------------------------------------------------
                    cond6 = 0;%ampR <= 0.5; % Amplitude ration condition
                    if cond6
                        continue;
                    end
                    %------------------------------------------------------
                    cond4 = (m_cr >= crESW) || (m_crdl >= crEDL); % condition for topology ESW or EDL
                    if cond4
                        %------------------------------------------------------
                        cond5 = T0 > Rbold ; % condition to avoid double counting
                        if cond5
                            %------------------------------------------------------
                            plot(ts(1:end-1)+lagsT(ncr)*dt-ts(1), Et', '--b', 'LineWidth', lw)
                            hold on
                            plot(ts(1:end-1)-ts(1), Es, 'r', 'LineWidth', lw)
                            
                            %                         plot(dx, Es, 'linewidth',lw)
                            
                            
                            
                            plot([ts(lb)+lagsT(ncr)*dt, ts(lb)+lagsT(ncr)*dt]-ts(1), [-1*ampH, ampH],'--k', 'LineWidth', lw-1)
                            plot([ts(rb)+lagsT(ncr)*dt, ts(rb)+lagsT(ncr)*dt]-ts(1), [-1*ampH, ampH],'--k', 'LineWidth', lw-1)
                            plot([ts(t0c)+lagsT(ncr)*dt, ts(t0c)+lagsT(ncr)*dt]-ts(1), [-1*ampH, ampH],'--m', 'LineWidth', lw-1)
                            ylabel(gca,'$E\ (mV/m)$','interpreter','latex','fontsize',15)
                            s=seconds(ts(1));
                            s.Format = 'hh:mm:ss.SSS';
                            xtring = ['Time after ', char(s)];
                            xlabel(gca,xtring,'interpreter','latex','fontsize',15)
                            
                            %irf_legend(gca,'(b)',[0.99 0.98],'color','k','fontsize',12)
                            
                            %tstring = irf_time(tm(J), 'epoch>utc');
                            %title(['t0 = ', tstring], 'FontSize',fsiz)
                            grid on
                            ax = gca;
                            ax.GridLineStyle = '--';
                            set(gca,'FontSize',fsiz,'FontWeight','bold','linewidth',lw-1)
                            
                            %
                            fname = [figpath, 'ESW_',sprintf('%05d',KJ)];
                            figsave(fname, 300, [8 6], 'png')
                            pause(0.02)
                            clf
                            %------------------------------------------------------
                            %                         pause
                            %                         close all;
                            SORTM(KJ, 1) = tm(Lb);
                            SORTM(KJ, 2) = tm(Rb);
                            SORTM(KJ, 3) = tm(T0);
                            SORTM(KJ, 4) = wd;
                            KJ = KJ+1;
                            %fprintf([repmat('\b', 1, 10) '%4.0f\n'],KJ);
                            
                        end
                    end
                    Rbold = Rb;
                end
            end
        end
    end
end

end