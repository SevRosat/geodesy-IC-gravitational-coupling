clear all;clc;format long;
Cm = 7.129*1e37; % kg.m2
icorr = 1; % = 1 if we correct from hydro (using ISBA) S22, else icorr = 0
% LOD is corrected from AOHSl
%==========================================================================
%               BAYESIAN INVERSION OF THE IC DIFFERENTIAL ROTATION
%                                   S. Rosat
%                                               May 2026
%==========================================================================
% Input: S22, alpha ~ dLOD/dt
% Parameters to inverse:
% gamma (gravi. torque), dh (topo ICB)
% Time unit: year
%==========================================================================
% Data
% S22 from SLR corrected for ISBA
% months since 2002
cfile = 'cs22_ts_new.nc';
% time_grace = ncread(cfile,'time_grace'); % 'hours since 2002-01-16 12:00:00'
% time_monthly  = ncread(cfile,'time_igg'); % 'hours since 1984-01-16 12:00:00'
% time_isba  = ncread(cfile,'time_isba'); % 'hours since 1979-01-16 12:00:00'
time_monthly  = ncread(cfile,'time_monthly'); % 'hours since 1984-01-16 12:00:00'
% s22_n = ncread(cfile,'s22_filt_igg_slr_doris');
% s22_n = ncread(cfile,'s22_igg_slr_upwr');
% s22_n = ncread(cfile,'s22_filt_igg_slr_doris');
s22_n = ncread(cfile,'s22_filt_igg_slr_hydrid');
% s22_n = ncread(cfile,'s22_filt_grace_csr');
if icorr==1
s22_isba = ncread(cfile,'s22_filt_isba');
s22_n = s22_n-s22_isba;
end
time_monthly_mat = datetime(time_monthly*3600,'ConvertFrom','epochtime','Epoch','1984-01-16 12:00:00');
% remove NaN
s22 = s22_n(~isnan(s22_n));
time_monthly_mat = time_monthly_mat(~isnan(s22_n));
%
% Take a shorter time window for tests
% s22 = s22(100:end);
% time_monthly_mat = time_monthly_mat(100:end);
%
sigma22 = std(s22);
%
% dLOD from IERS excitation functions corrected for AOHSl
% The data is alpha with gamma = 1
% Open the LOD file in ms.
% filename = 'B:\TRAVAUX_RECHERCHE\Differential_IC_rotation\Chiz1976-2026_AOHSlcorr.txt';
% startRow = 4;
% formatSpec = '%9f%9f%10f%10f%10f%[^\n\r]';
% fileID = fopen(filename,'r');
% textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
% dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% for block=2:length(startRow)
%     frewind(fileID);
%     textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
%     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%     for col=1:length(dataArray)
%         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%     end
% end
% fclose(fileID);
% Chiz19762026 = table(dataArray{1:end-1}, 'VariableNames', {'yrs','mjd','chix','chiy','chiz'});
% % Conversion des années en dates MATLAB (pour faciliter les calculs de différences)
% dateschi = decyr2datenum(Chiz19762026.yrs);
% % filtering
fs = 1/0.0027; % in cpy
% Tfi = 3; % in years
% norder = 3;
% % LP filtering
% [b,a] = butter(norder,1/Tfi/fs*2, 'low');
% Chizfi = filtfilt(b,a,Chiz19762026.chiz);
% % sigmaLOD = std(Chizfi)./1000/86400;
%
% LOD deconvolved from Core 
% load('B:\TRAVAUX_RECHERCHE\ROTATION\LOD\Torques\LODdeconv-GGJx22.mat');
% ATTENTION sampling rate different --> to be modified
load('B:\TRAVAUX_RECHERCHE\ROTATION\LOD\Torques\LODdeconv-IGFx23.mat');
%'newYears','LODres'
% resample to fs = 1 month
dateschi0 = decyr2datenum(newYears);
[Chizfi,dateschi] = resample(LODres,dateschi0,fs);
%
% alpha angle
alpha1 = -2*pi/86400^2*Cm*diff(Chizfi/1000)./diff(dateschi*86400);
% smoothing
% alpha1 = sgolayfilt(alpha1,3,501);
talpha1 = datetime(dateschi(1:end-1),'ConvertFrom','Datenum');
% interpolate alpha1 at S22 times
alpha_timeS22 = interp1(talpha1, alpha1, time_monthly_mat,'linear','extrap');
%
% figure
% plot(time_monthly_mat,alpha1_new,'r');
% hold on
% plot(talpha1,alpha1,'k');
%
% sigma = sqrt(sigma22^2+sigmaLOD^2);
sigma = sigma22;
% Number of points for one side of the parameter space
%==========================================================================
nb = 50;
nc = 50;
nd = 50;
%
% Parameters space
%==========================================================================
dh = linspace(10,100,nb); % m
% x = log10(gamma)
x = linspace(18,23,nc); 
gamma = 10.^x; % N.m
tau = linspace(0.01,10,nd); % year
%
% matrix initialization
%==========================================================================
rhoD = zeros(nb,nc,nd);
% margdhG = zeros(nb,nc);
% margdhtau = zeros(nb,nd);
% margGtau = zeros(nc,nd);
%
% a priori probability density
%==========================================================================
rhoM = 1;
%
% Waitbar initialization
%==========================================================================
tic
h = waitbar(0,'Processing Bayesian inversion, please wait...');
%
% Computation of the G(m) induced by the IC reorientation
% and computation of the a posteriori probability density
%==========================================================================
iw = 0;
for ii = 1:nb % dh
    for jj = 1:nc % gamma
        for kk = 1:nd % tau
        % phi angle
        phim = alpha1./gamma(jj)-2*pi/86400^2*Cm./gamma(jj)*Chizfi(1:end-1)/1000./tau(kk)/365.25/86400;
        if length(phim)>1000 && ii==1 && jj==1 && kk==1
            disp('smoothing');
            phim = sgolayfilt(phim,3,501);
        end
        phi_timeS22 = interp1(talpha1, phim, time_monthly_mat,'linear','extrap');
        alpha = phi_timeS22 + (alpha_timeS22(1)-phi_timeS22(1)).*exp(-(datenum(time_monthly_mat))./tau(kk)/365.25);
        s22m = 1e-11*dh(ii)*alpha;
        % Resonance function
        rhoD(ii,jj,kk) = exp(-0.5*sum(abs(s22m-s22).^2./abs(sigma).^2));
        iw = iw + 1;
        % Waitbar update
        waitbar(iw/nb/nc/nd,h);
        end
    end
end
%
% a posteriori pdf
pdf = rhoM.*rhoD;
%
% Computation of marginal laws
%==========================================================================
% dh, gamma, tau
margdhG = sum(pdf,3);
margdhG = margdhG./norm(margdhG);
margGtau = squeeze(sum(pdf,1));
margGtau = margGtau./norm(margGtau);
margdhtau = squeeze(sum(pdf,2));
margdhtau = margdhtau./norm(margdhtau);
% marginal laws
margG2 = sum(margdhG,1);
margG = sum(margGtau,2);
margtau = sum(margGtau,1);
margtau2 = sum(margdhtau,1);
%
% Close waitbar
%==========================================================================
close(h);
% elapsed time
toc
%
% Plots 
%==========================================================================
% x-T
subplot(1,3,1)
imagesc(dh,x,margdhG.');axis xy;
xlabel('$\delta$h (m)','fontsize',16);
ylabel('log($\Gamma$) (N.m)','fontsize',16);
set(gca,'fontsize',14);
axis('square');
title('(a)')

subplot(1,3,2)
imagesc(dh,tau,margdhtau.');axis xy;
xlabel('$\delta$h (m)','fontsize',16);
ylabel('$\tau$ (year)','fontsize',16);
set(gca,'fontsize',14);
axis('square');
title('(b)')
% 
subplot(1,3,3);
imagesc(x,tau,margGtau.');axis xy;
xlabel('log($\Gamma$) (N.m)','fontsize',16);
ylabel('$\tau$ (year)','fontsize',16);
set(gca,'fontsize',14);
axis('square');
title('(c)')
%
% normdh = sum(margdh);
normG = sum(margG);
normtau = norm(margtau);
% dhmean = sum(margdh/normdh.*dh.');
xmean = sum(margG.'/normG.*x);
Gmean = 10^xmean;
taumean = sum(margtau/normtau);
% vardh = sum(margdh/normdh.*(dh.'-dhmean).^2);
varx = sum(margG.'/normG.*(x-xmean).^2);
% varG = 10^varx;
vartau = sum(margtau/normtau.*(tau-taumean).^2);

disp('');
disp(['Mean value of tau = ',num2str(taumean), ' year', ' std = ', num2str(sqrt(vartau))]);
%num2str(10^Xmean*log(10)*sqrt(varX))
disp(['Mean value of Gamma = ',num2str(Gmean),' N.m',' std = ', num2str(10^xmean*log(10)*sqrt(varx))]);
disp('');

figure
subplot(1,2,1)
% plot(dh,margdh.','k');
% xlabel ('\deltah (m)','fontsize',20)
% set(gca,'fontsize',14);

semilogx(gamma,margG,'k');
hold on
% semilogx(gamma,margG2);
plot([Gmean Gmean],[0 max(margG)],'k--')
xlabel ('$\Gamma$ (N.m)','fontsize',20)
set(gca,'fontsize',14);
axis 'square'
title('(a)')
%
subplot(1,2,2)
plot(tau,margtau.','k');
hold on
plot([taumean taumean],[0 max(margtau)],'k--')
xlabel ('$\tau$ (year)','fontsize',20)
set(gca,'fontsize',14);
axis 'square'
title('(b)')