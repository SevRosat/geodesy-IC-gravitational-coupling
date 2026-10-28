clear; clc; format long;
%==========================================================================
%               BAYESIAN INVERSION OF THE IC DIFFERENTIAL ROTATION
%                                   S. Rosat
%                     May 2026 - corrected version (oct. 2026)
%                     + 1-sigma uncertainties on the MAP solution
%==========================================================================
% Parameters : dh (topo ICB), Gamma (gravitational coupling), tau (visco)
% ATTENTION : S22 model depends only on dh/Gamma and tau,
% dh and Gamma are then non separable without independent a priori.
%==========================================================================
Cm    = 7.129e37;   % kg.m2
icorr = 0;          % 1 : S22 corrected from ISBA hydrology
sigma_n = [];       % standard deviation of SLR series (ierr = 0 only)
ierr  = 2;          % data errors in the likelihood:
                    %   0 : constant sigma (sigma_n, or std(s22) if empty)
                    %   1 : sigma_i = e_i, error bars of the series taken at face value
                    %   2 : sigma_i = s*e_i, e_i give the relative weights and the
                    %       scale factor s is unknown and marginalized analytically
                    %       (Jeffreys prior 1/s) -> robust if errors are under/over-estimated

%% S22 ---------------------------------------------------------------------
cfile = 'cs22_ts_errors.nc';
sname = 's22_filt_igg_slr_hydrid';
% sname = 's22_igg_slr_upwr';
% sname = 's22_filt_igg_slr_doris';
% sname = 's22_filt_grace_csr';
s22_n = double(ncread(cfile,sname));  s22_n = s22_n(:);
% error bars of the series (variable 'e' + name in the netCDF file)
info_nc = ncinfo(cfile);
if any(strcmp({info_nc.Variables.Name},['e' sname]))
    e22_n = double(ncread(cfile,['e' sname]));  e22_n = e22_n(:);
elseif ierr > 0
    error('No error variable e%s in %s: use ierr = 0',sname,cfile);
else
    e22_n = nan(size(s22_n));
end
if icorr==1
    % ISBA has no error estimate: the error of the corrected series is that of SLR
    s22_n = s22_n - double(reshape(ncread(cfile,'s22_filt_isba'),[],1));
end
% time axis from its 'units' attribute
% (in cs22_ts_errors.nc: 'days since 1984-01-01')
t22 = nc_time(cfile,'time_monthly');  t22 = t22(:);
ok  = ~isnan(s22_n);
if ierr > 0, ok = ok & ~isnan(e22_n) & e22_n > 0; end
s22 = s22_n(ok);  e22 = e22_n(ok);  t22 = t22(ok);
%
N   = numel(s22);
if isempty(sigma_n), sigma = std(s22); else, sigma = sigma_n; end
% weights w_i = 1/sigma_i^2 (relative weights only for ierr = 2)
switch ierr
    case 0, w = ones(N,1)/sigma^2;
    case {1,2}, w = 1./e22.^2;
    otherwise, error('ierr must be 0, 1 or 2');
end
fprintf('N = %d data, mean error = %.2e, std(S22) = %.2e\n', N, sqrt(mean(e22.^2)), std(s22));

%% LOD (excitation chi_z corrected from AOHSl, in ms) ------------------------------
filename = 'Chiz1976-2026_AOHSlcorr.txt';              
D    = readmatrix(filename,'FileType','text','NumHeaderLines',3);
chiz = D(:,5);                                         % ms
mjd  = D(:,2);
tL   = datetime(mjd,'ConvertFrom','mjd');              % MJD  
                                                      
dt   = 86400;                                          
fs   = 365.25;                                         
Tfi  = 3;  norder = 3;
[b,a] = butter(norder, (1/Tfi)/(fs/2), 'low');
Lf   = filtfilt(b,a,chiz)/1000;                        % s

K      = 2*pi/86400^2*Cm;
KL     = -K*Lf;   
torque = -K*gradient(Lf, dt);                          % N.m

A = interp1(tL, torque, t22);
B = interp1(tL, KL,     t22);                        % N.m (to be divided by tau in years)
tyr = days(t22 - t22(1))/365.25;                       

%% Parameter space -----------------------------------------------------
nb = 100; nc = 100; nd = 100;
dh    = linspace(10,100,nb);       % m
x     = linspace(18,23,nc);        % log10(Gamma)
gamma = 10.^x;                     % N.m
tau   = linspace(0.1,20,nd);      % year

%% Log-likelihood -----------------------------------------------------------
chi2 = inf(nb,nc,nd);                % weighted misfit sum(w_i*r_i^2)
tic
for kk = 1:nd
    E = exp(-tyr/tau(kk));
    for jj = 1:nc
        phi   = (A + B/tau(kk)./365.25/86400)/gamma(jj);        % rad
        alpha = phi + (A(1)/gamma(jj) - phi(1)).*E;             % alpha
        r     = s22 - 1e-11*alpha*dh;                           % N x nb
        chi2(:,jj,kk) = sum(w.*r.^2,1).';
    end
end
toc
if ierr < 2
    logL = -0.5*fcorr*chi2;
else
    % int s^-Neff exp(-fcorr*chi2/(2 s^2)) ds/s  propto  chi2^(-Neff/2), Neff = fcorr*N
    logL = -0.5*fcorr*N*log(chi2);
end

pdf = exp(logL - max(logL(:)));                                 % no underflow
pdf = pdf/sum(pdf(:));

%% Marginal laws (normalized by the SUM) -------------------------------------
margdhG   = sum(pdf,3);                    % nb x nc
margdhtau = squeeze(sum(pdf,2));           % nb x nd
margGtau  = squeeze(sum(pdf,1));           % nc x nd
margdh  = sum(margdhG,2);                  % nb x 1
margG   = sum(margdhG,1).';                % nc x 1
margtau = sum(margdhtau,1);                % 1 x nd

dhmean  = sum(margdh.*dh(:));    vardh  = sum(margdh.*(dh(:)-dhmean).^2);
xmean   = sum(margG.*x(:));      varx   = sum(margG.*(x(:)-xmean).^2);
taumean = sum(margtau.*tau);     vartau = sum(margtau.*(tau-taumean).^2);   

% The only quantity that is really resolved by the data: q = log10(dh/Gamma)
Q     = log10(dh(:)) - x;                  % nb x nc
qmean = sum(margdhG.*Q,'all');  qstd = sqrt(sum(margdhG.*(Q-qmean).^2,'all'));

chi2min = min(chi2(:));
s22w    = sum(w.*s22)/sum(w);                    % weighted mean
chi2_0  = sum(w.*(s22-s22w).^2);
Neff    = fcorr*N;
fprintf('\nVariance reduction at MAP : %.1f %%\n', 100*(1-chi2min/chi2_0));
if ierr == 1
    % ~1 if the error bars are consistent with the misfit; >>1 : errors underestimated
    fprintf('Reduced chi2 at MAP : %.2f  (>>1 : error bars too small, try ierr = 2)\n', ...
        fcorr*chi2min/(Neff-3));
elseif ierr == 2
    fprintf('Error scale factor at MAP : s = %.2f\n', sqrt(chi2min/(N-3)));
end
fprintf('tau   = %.2f +/- %.2f year\n', taumean, sqrt(vartau));
fprintf('log10(Gamma) = %.2f +/- %.2f  (Gamma ~ %.2e ×/÷ %.0f N.m)\n', xmean, sqrt(varx), 10^xmean, 10^sqrt(varx));
fprintf('dh    = %.1f +/- %.1f m  (not constrained alone)\n', dhmean, sqrt(vardh));
fprintf('log10(dh/Gamma) = %.2f +/- %.2f   <- resolved parameter\n\n', qmean, qstd);

%% MAP (maximum of the joint posterior) ---------------------------------------
% With flat priors on (dh, log10Gamma, tau) the MAP is the max of logL.
[~,imap]       = max(logL(:));
[ib,ic,id]     = ind2sub(size(logL),imap);
dh_map    = dh(ib);
x_map     = x(ic);       gamma_map = gamma(ic);
tau_map   = tau(id);
q_map     = log10(dh_map) - x_map;
chi2_map  = chi2(ib,ic,id);                      % = chi2min
onedge    = [ib==1||ib==nb, ic==1||ic==nc, id==1||id==nd];

% Maxima of the 1D marginals (can differ from the joint MAP)
[~,i1] = max(margdh);   dh_mode  = dh(i1);
[~,i2] = max(margG);    x_mode   = x(i2);
[~,i3] = max(margtau);  tau_mode = tau(i3);

% Model at the MAP
E_map     = exp(-tyr/tau_map);
phi_map   = (A + B/tau_map/365.25/86400)/gamma_map;            % rad
alpha_map = phi_map + (A(1)/gamma_map - phi_map(1)).*E_map;    % rad
omega_map = gradient(phi_map,tyr);                             % rad/yr
s22_map   = 1e-11*alpha_map*dh_map;
res_map   = s22 - s22_map;

% Error bars used in the likelihood, at the MAP
switch ierr
    case 0, sig_plot = sigma*ones(N,1);
    case 1, sig_plot = e22;
    case 2, s_map = sqrt(chi2_map/(N-3)); sig_plot = s_map*e22;
end

%% 1-sigma uncertainties on the MAP (profile likelihood) ----------------------
% 1-sigma interval of a parameter = range where its profile likelihood
% (logL maximized over the other parameters) stays within 0.5 of the maximum.
Lp       = max(logL,[],3);                           % nb x nc (max over tau)
prof_dh  = max(Lp,[],2);                             % nb x 1
prof_x   = max(Lp,[],1).';                           % nc x 1
prof_tau = squeeze(max(max(logL,[],1),[],2));        % nd x 1
% profile of q (logL depends on dh and Gamma only through q)
qg     = linspace(min(Q(:)), max(Q(:)), 300);
prof_q = -inf(size(qg));
for i = 1:nb
    prof_q = max(prof_q, interp1(x, Lp(i,:), log10(dh(i)) - qg, 'linear', -inf));
end
[dh_lo,  dh_hi,  e_dh ] = sig1(dh,  prof_dh);
[x_lo,   x_hi,   e_x  ] = sig1(x,   prof_x);
[tau_lo, tau_hi, e_tau] = sig1(tau, prof_tau);
[q_lo,   q_hi,   e_q  ] = sig1(qg,  prof_q);
edgetxt = {'',' (bounded by the grid only)'};

%% Print MAP +/- 1 sigma ------------------------------------------------------
fprintf('------------------ MAP (+/- 1 sigma) ------------------\n');
fprintf('dh_map           = %.1f  (+%.1f / -%.1f) m%s\n', dh_map, dh_hi-dh_map, dh_map-dh_lo, edgetxt{e_dh+1});
fprintf('log10(Gamma_map) = %.2f (+%.2f / -%.2f)%s\n', x_map, x_hi-x_map, x_map-x_lo, edgetxt{e_x+1});
fprintf('                   Gamma_map = %.2e N.m  [%.2e, %.2e]\n', gamma_map, 10^x_lo, 10^x_hi);
fprintf('tau_map          = %.2f (+%.2f / -%.2f) year%s\n', tau_map, tau_hi-tau_map, tau_map-tau_lo, edgetxt{e_tau+1});
fprintf('log10(dh/Gamma)_map = %.2f (+%.2f / -%.2f)%s   <- resolved parameter\n', q_map, q_hi-q_map, q_map-q_lo, edgetxt{e_q+1});
fprintf('chi2 at MAP      = %.4e   (wrms residual = %.2e, rms residual = %.2e)\n', ...
    chi2_map, sqrt(chi2_map/sum(w)), sqrt(mean(res_map.^2)));
fprintf('Modes of 1D marginals: dh = %.1f m, log10(Gamma) = %.2f, tau = %.2f yr\n', ...
    dh_mode, x_mode, tau_mode);
names = {'dh','log10(Gamma)','tau'};
for p = find(onedge)
    fprintf('WARNING: MAP of %s is on the edge of the grid -> extend the range\n', names{p});
end
fprintf('-------------------------------------------------------\n\n');

%% Figures --------------------------------------------------------------------
set(groot,'defaultTextInterpreter','latex');            
set(groot,'defaultAxesTickLabelInterpreter','latex');

% 2D marginals, with joint MAP (white +) and posterior means (white o)
figure
subplot(1,3,1); imagesc(dh,x,margdhG.'); axis xy square; hold on
plot(dh, log10(dh) - q_map, 'k--', 'LineWidth', 1)          % ridge dh/Gamma = MAP value
plot(dh, log10(dh) - q_lo, 'k:', 'LineWidth', 0.8)   % 1-sigma bounds on q
plot(dh, log10(dh) - q_hi, 'k:', 'LineWidth', 0.8)
plot(dh_map, x_map, 'w+', 'MarkerSize', 10, 'LineWidth', 1.5)
plot(dhmean, xmean, 'wo', 'MarkerSize', 6, 'LineWidth', 1)
ylim([x(1) x(end)])
xlabel('$\delta h$ (m)'); ylabel('$\log_{10}\Gamma$ (N.m)'); title('(a)')
subplot(1,3,2); imagesc(dh,tau,margdhtau.'); axis xy square; hold on
plot(dh_map, tau_map, 'w+', 'MarkerSize', 10, 'LineWidth', 1.5)
plot(dhmean, taumean, 'wo', 'MarkerSize', 6, 'LineWidth', 1)
xlabel('$\delta h$ (m)'); ylabel('$\tau$ (yr)'); title('(b)')
subplot(1,3,3); imagesc(x,tau,margGtau.'); axis xy square; hold on
plot(x_map, tau_map, 'w+', 'MarkerSize', 10, 'LineWidth', 1.5)
plot(xmean, taumean, 'wo', 'MarkerSize', 6, 'LineWidth', 1)
xlabel('$\log_{10}\Gamma$ (N.m)'); ylabel('$\tau$ (yr)'); title('(c)')

% 1D marginals: mean (dashed), joint MAP (red), MAP 1-sigma interval (red dotted)
figure
subplot(1,2,1)
semilogx(gamma,margG,'k'); hold on
plot(10^xmean*[1 1],[0 max(margG)],'k--')
plot(gamma_map*[1 1],[0 max(margG)],'r-')
plot(10.^[x_lo x_lo NaN x_hi x_hi],[0 max(margG) NaN 0 max(margG)],'r:')
legend('marginal','mean','MAP','MAP $\pm1\sigma$','Interpreter','latex','Location','best')
xlabel('$\Gamma$ (N.m)'); axis square; title('(a)')
subplot(1,2,2)
plot(tau,margtau,'k'); hold on
plot(taumean*[1 1],[0 max(margtau)],'k--')
plot(tau_map*[1 1],[0 max(margtau)],'r-')
plot([tau_lo tau_lo NaN tau_hi tau_hi],[0 max(margtau) NaN 0 max(margtau)],'r:')
legend('marginal','mean','MAP','MAP $\pm1\sigma$','Interpreter','latex','Location','best')
xlabel('$\tau$ (yr)'); axis square; title('(b)')

% Data fit at the MAP
figure
subplot(3,1,1:2)
errorbar(t22, s22, sig_plot, '.', 'Color', [0.6 0.6 0.6], 'CapSize', 0); hold on
plot(t22, s22, 'k.')
plot(t22, s22_map, 'r', 'LineWidth', 1.5)
ylabel('$S_{22}$'); grid on; title('(a)')
legend('data $\pm\sigma$','data', ...
    sprintf('MAP: $\\tau$=%.1f$^{+%.1f}_{-%.1f}$ yr, $\\log_{10}(\\delta h/\\Gamma)$=%.2f$^{+%.2f}_{-%.2f}$', ...
    tau_map, tau_hi-tau_map, tau_map-tau_lo, q_map, q_hi-q_map, q_map-q_lo), ...
    'Interpreter','latex','Location','best')
subplot(3,1,3)
plot(t22, res_map, 'k.-'); hold on
plot(t22, sig_plot, 'r:', t22, -sig_plot, 'r:')
ylabel('residual'); grid on; title('(b)')

% Yang & Song (2023, Nat. Geosci., doi:10.1038/s41561-022-01112-z) -----------
% Schematic sinusoidal version of their multidecadal oscillation:
% differential rotation rate vanishing at the turning
% points of the early 1970s and ~2009 (pause), period P = 2*(t2-t1) ~ 7 decades,
% amplitude set so that the mean rate over 1991-2010 is 0.127 deg/yr
% (Yang & Song 2022, GRL), angle set to 0 deg in 2000 (as in their Fig. 4).
ys.t1 = 1972;  ys.t2 = 2009;  ys.P = 2*(ys.t2 - ys.t1);   % yr
ys.Om_ref = 0.127;  ys.win = [1991 2010];  ys.tz = 2000;  % deg/yr, yr, yr
f_ys  = @(t) sin(2*pi*(t - ys.t1)/ys.P);                   % normalized rate
tw_ys = linspace(ys.win(1), ys.win(2), 1000);
ys.Oma = ys.Om_ref/mean(f_ys(tw_ys));                      % rate amplitude (deg/yr)
Om_ys  = @(t) ys.Oma*f_ys(t);                              % deg/yr
phi_ys = @(t) -ys.Oma*ys.P/(2*pi)*(cos(2*pi*(t - ys.t1)/ys.P) - cos(2*pi*(ys.tz - ys.t1)/ys.P)); % deg
fprintf('Yang & Song (2023), schematic: rate amplitude %.3f deg/yr, period %d yr\n', ys.Oma, ys.P);
% comparison on Omega only (the angle phi_ys is defined up to an arbitrary offset)
ys_omega = Om_ys(tyr);
% Phi-angle from Yang&Song 2023
temp = load("ICangle_from_seismology1995-2020_YangSong2023.txt");
tsismo = datetime(datenum(temp(:,1),0,0),'ConvertFrom','datenum');
phisismo = temp(:,2); % in degree

% Rotation angle of the IC at the MAP
figure
subplot(3,1,1)
plot(t22, phi_map*180/pi, 'k', 'LineWidth', 1.2); grid on
hold on
plot(tsismo, phisismo,'r-', 'linewidth', 2);
ylabel('$\phi$ ($^\circ$)')
title(sprintf('IC differential rotation at MAP ($\\tau$ = %.1f yr, $\\log_{10}(\\delta h/\\Gamma)$ = %.2f)', tau_map, q_map))

subplot(3,1,2)
plot(t22, alpha_map*180/pi, 'k', 'LineWidth', 1.2); grid on
ylabel('$\alpha$ ($^\circ$)')

subplot(3,1,3)
plot(t22,omega_map*180/pi,'k', 'LineWidth', 1.2); grid on
hold on
plot(t22, ys_omega, 'r-', 'linewidth', 2);
xlabel('Time'); ylabel('$\Omega$ ($^\circ$/yr)')

%% Local functions ------------------------------------------------------------
function t = nc_time(cfile,var)
% netCDF time variable -> datetime, using its 'units' attribute
% ('<days|hours|minutes|seconds> since yyyy-MM-dd [HH:mm:ss]')
t = double(ncread(cfile,var));
u = strtrim(ncreadatt(cfile,var,'units'));
parts = strsplit(u);
unit  = lower(parts{1});
ref   = strjoin(parts(3:end),' ');
fmts  = {'yyyy-MM-dd HH:mm:ss','yyyy-MM-dd HH:mm:ss.S','yyyy-MM-dd HH:mm','yyyy-MM-dd'};
t0 = NaT;
for f = 1:numel(fmts)
    try
        t0 = datetime(ref,'InputFormat',fmts{f});
        break
    catch
    end
end
if isnat(t0), error('Cannot parse time reference "%s" of %s',ref,var); end
switch unit
    case {'days','day'},       t = t0 + days(t);
    case {'hours','hour'},     t = t0 + hours(t);
    case {'minutes','minute'}, t = t0 + minutes(t);
    case {'seconds','second'}, t = t0 + seconds(t);
    otherwise, error('Unknown time unit "%s" for %s',unit,var);
end
end

function [lo,hi,edge] = sig1(xg, Lp)
% 1-sigma profile-likelihood interval {x : logL_prof(x) >= max - 0.5},
% on a 20x finer grid; edge = true if the interval reaches the grid boundary
xf = linspace(xg(1), xg(end), 20*numel(xg));
Lf = interp1(xg(:), Lp(:), xf, 'linear');
ok = Lf >= max(Lf) - 0.5;
lo = min(xf(ok));  hi = max(xf(ok));
edge = ok(1) || ok(end);
end
