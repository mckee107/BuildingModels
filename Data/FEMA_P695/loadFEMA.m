% Codi McKee
% Texas A&M University
% 02-Jul-2020

% FEMA P695 Far-Field Ground Motions

clear all; clc; close all
addpath([pwd '\..\..\Other Functions'])

% Load Data
load FEMA_P695.mat;
load TargetSpectra.mat

if exist('GroundMotions.mat', 'file') == 2
    load GroundMotions.mat
else
    GroundMotions = struct('No',[],'Event',[],'Date',[],'Station',[],'Component',[],'Units',[],'NPTS',[],'DT',[],'TimeSeries',[]);
    save GroundMotions.mat GroundMotions
end

% Plot Scaled Spectra

figure(1); hold all;
for ii = 1:22
   plot(FEMA_P695.FarField(ii).Scaled_Spectra_TXYZ(:,1),FEMA_P695.FarField(ii).Scaled_Spectra_TXYZ(:,2),'Color',[0.50 0.50 0.50])
end
% set(gca,'Xscale','log')
xlim([0 4])

% Geo Mean
geoMean = zeros(length(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1)),1);
mMean = zeros(length(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1)),1);


for ii = 1:length(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1))
    geoMean(ii) = 1;
    
    for jj = 1:22
        for kk = 1:2
            geoMean(ii) = geoMean(ii) * FEMA_P695.FarField(jj).Scaled_Spectra_TXYZ(ii,kk+1);
            mMean(ii)   = mMean(ii) + FEMA_P695.FarField(jj).Scaled_Spectra_TXYZ(ii,kk+1);
        end
    end
    geoMean(ii) = geoMean(ii) ^ (1/44);
    mMean(ii) = mMean(ii)/44;
    
end
plot(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1),geoMean,'Color',[0.0 0.0 0.0])
% plot(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1),mMean,'Color','r')
plot(TargetSpectra(:,1),TargetSpectra(:,2),'k')

% Scale ground motion suite

Tn = 0.4585738086971623;

geoMean_wn = interp1(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1),geoMean,Tn);

ts_wn = interp1(TargetSpectra(:,1),TargetSpectra(:,2),Tn);

scaleFac = ts_wn/geoMean_wn;

geoMean_new = geoMean .* scaleFac;
% plot(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1),geoMean_new,'Color','r')

figure(2); FigSize(4,4); hold all
plot([NaN NaN],[NaN NaN],'Color',[0.75 0.75 0.75],'LineWidth',0.5)
plot([NaN NaN],[NaN NaN],'Color','r','LineWidth',1.0)
plot([NaN NaN],[NaN NaN],'Color','k','LineWidth',1.0)
plot([NaN NaN],[NaN NaN],'--r','LineWidth',1.0)
plot([NaN NaN],[NaN NaN],'--k','LineWidth',1.0)
plot([NaN NaN],[NaN NaN],'Color','b','LineWidth',1.0)

for jj = 1:22
    for kk = 1:2
        GM = (jj-1)*2 + kk;
        GroundMotions(GM).DT = FEMA_P695.FarField(jj).dt;  
        GroundMotions(GM).NPTS = length(FEMA_P695.FarField(jj).Scaled_Motion_XYZ(:,1));  
        GroundMotions(GM).No = FEMA_P695.FarField(jj).ID; 
        GroundMotions(GM).TimeSeries = FEMA_P695.FarField(jj).Scaled_Motion_XYZ(:,kk) .* scaleFac; 
        GroundMotions(GM).Spectra(:,1) = FEMA_P695.FarField(jj).Scaled_Spectra_TXYZ(:,1);
        GroundMotions(GM).Spectra(:,2) = FEMA_P695.FarField(jj).Scaled_Spectra_TXYZ(:,kk+1).* scaleFac;
        
        plot(GroundMotions(GM).Spectra(:,1),GroundMotions(GM).Spectra(:,2),'Color',[0.75 0.75 0.75],'LineWidth',0.5)
    end
end
plot(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1),geoMean_new,'Color','r','LineWidth',1.0)
plot(TargetSpectra(:,1),TargetSpectra(:,2),'k','LineWidth',1.0)
plot(FEMA_P695.FarField(1).Scaled_Spectra_TXYZ(:,1),geoMean_new.*1.5,'--r','LineWidth',1.0)
plot(TargetSpectra(:,1),TargetSpectra(:,2).*1.5,'--k','LineWidth',1.0)
xlim([0 4])
yLims = ylim;
plot([Tn Tn],[yLims(1) yLims(2)],'-b','LineWidth',1.0)
legend('Scaled Spectra','DE Geometric Mean','DE Target Spectra','MCE Geometric Mean','MCE Target Spectra','T_N_1')
title('Scaled FEMA P695 Spectra')
xlabel('Period [s]')
ylabel('S_a [g]')
FigureFormatFunction

tightfig; FigSize(4,4)

% if printVal
%     export_fig([pwd '\Scaled Spectra'], '-png', '-r1200','-nocrop')
%     savefig([pwd '\Scaled Spectra'])
% end

figure; hold all
for GM = 1:2

t = [0:GroundMotions(GM).DT:(GroundMotions(GM).NPTS-1)*GroundMotions(GM).DT]; 
plot(t,GroundMotions(GM).TimeSeries);

end

save GroundMotions.mat GroundMotions 


