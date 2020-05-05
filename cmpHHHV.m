%% Common Midpoint Gather Dual-Pole (HH-HV) GPR (cmpHHHV)
clear; close all; clc;
%% Meta Data is a structure MD
MD.dataDir = 'Enter Data Directory Here.';

% Import Utilities
addpath './functions';
addpath './colormaps';
isSave = 0;

% Extract Files
MD.workDir = pwd;
MD.fileNames = dir([MD.dataDir,'/','*.nc']);
MD.lineNo = 1:length(MD.fileNames);   % Array of data "LINE" numbers
MD.nFiles = length(MD.lineNo);        % Number of Files

%% Processing WorkFlow Controls
% Read Data
isReadNC = 1;           % Read and Process CMP Data
isTrimTWT = 0;          % Truncate Recorded Data
isReduceData = 0;       % Remove Traces
isVelocityAnalysis = 1; % Picking and Inversion
%% Control Parallelization
% Parallel Computing Enable
isParallel = 1;
if isParallel
    % Wake Parallel Computing
    if isempty(gcp('nocreate'))     
        nWorkers = 2;
        p = parpool(nWorkers);
    else nWorkers = 2;
    end
else
    nWorkers = 1;
end

%% Read & Process GPR Data
if isReadNC
    processCMP_HHHV_NSIDC    
end
%% Velocity Analysis
if isVelocityAnalysis
    VelocityAnalysis
end
%% Save OutPut
if isSave
    save([fullfile(MD.fileNames(1).folder,MD.fileNames(1).name(1:end-6)),'.mat'],'MD','D')
end
%% Make Figures
    set(0,'DefaultAxesFontName','Serif')
    set(0,'DefaultTextFontName','Serif')
    
% Plot CMP Gather
rows = 1:151; % Trim Travel-time Axis
% Make Figure
figure();
imagesc(MD.offsetArray,MD.TimeAxis(rows),D.CMP{1}(rows,:)); colormap(bone);

% Plot Picks and Snow Parameters if isVelocityAnalysis
if isVelocityAnalysis
hold on;
for hh = 1:MD.ndir
plot(MD.offsetArray,D.dirTpick{hh},'rx','markersize',8,'linewidth',2)
end
for hh = 1:MD.nref
plot(MD.offsetArray,D.refTpick{hh},'rx','markersize',8,'linewidth',2)
end
titlestr = cell2mat(join(split(MD.fileNames(1).name(1:end-3),'_')));
title(titlestr)
xlabel('Offset (m)')
ylabel('Travel-time (ns)')
set(gca,'fontsize',14,'fontweight','bold')
daspect([0.3,1,1])

% Plot Parameter Distributions
% Direct Wave
subplotIx = 0;
for kk = 1:MD.ndir
subplotIx = subplotIx+1;
figure();
subplot(MD.ndir+MD.nref,2,subplotIx)
histogram(D.dirZ{kk},'facecolor','k','NumBins',10,'Normalization','probability')
title(['\fontsize{10}','t_0 = ',num2str(mean(D.dirT0{kk}),'%.1f'),' (ns)'])
ylabel('Probability')
xlabel('Depth (m)')
set(gca,'fontsize',10,'fontweight','bold','xtick',...
    [mean(D.dirZ{kk})-1.96.*std(D.dirZ{kk}),mean(D.dirZ{kk}),mean(D.dirZ{kk})+1.96.*std(D.dirZ{kk})])
xticklabels(round([mean(D.dirZ{kk})-1.96.*std(D.dirZ{kk}),mean(D.dirZ{kk}),mean(D.dirZ{kk})+1.96.*std(D.dirZ{kk})],2))
ylim([0 0.35])
sgtitle(titlestr,'fontweight','bold')
grid on

subplotIx = subplotIx+1;
subplot(MD.ndir+MD.nref,2,subplotIx)
histogram(D.dirRho{kk},'facecolor','k','NumBins',10,'Normalization','probability')
title({['\fontsize{10}','t_0 = ',num2str(mean(D.dirT0{kk}),'%.1f'),' (ns)']})
% ylabel('Probability')
xlabel('Density (g/cm^3)')
ylim([0 0.35])
% xlim([min(D.dirRho{kk}),max(D.dirRho{kk})])
set(gca,'fontsize',10,'fontweight','bold','xtick',...
    [mean(D.dirRho{kk})-1.96.*std(D.dirRho{kk}),mean(D.dirRho{kk}),mean(D.dirRho{kk})+1.96.*std(D.dirRho{kk})])
xticklabels(round([mean(D.dirRho{kk})-1.96.*std(D.dirRho{kk}),mean(D.dirRho{kk}),mean(D.dirRho{kk})+1.96.*std(D.dirRho{kk})],2))
grid on
end

% Reflections
for kk = 1:MD.nref
subplotIx = subplotIx+1;
subplot(MD.ndir+MD.nref,2,subplotIx)
histogram(D.refZ{kk},'facecolor','k','NumBins',10,'Normalization','probability')
xlabel('Depth (m)')
ylabel('Probability')
title({['\fontsize{10}','t_0 = ',num2str(mean(D.refT0{kk}),'%.1f'),' (ns)']})
set(gca,'fontsize',10,'fontweight','bold','xtick',...
    [mean(D.refZ{kk})-1.96.*std(D.refZ{kk}),mean(D.refZ{kk}),mean(D.refZ{kk})+1.96.*std(D.refZ{kk})])
xticklabels(round([mean(D.refZ{kk})-1.96.*std(D.refZ{kk}),mean(D.refZ{kk}),mean(D.refZ{kk})+1.96.*std(D.refZ{kk})],2))
ylim([0 0.35])
grid on

subplotIx = subplotIx+1;
subplot(MD.ndir+MD.nref,2,subplotIx)
histogram(D.refRho{kk},'facecolor','k','NumBins',10,'Normalization','probability')
xlabel('Density (g/cm^3)')
title(['\fontsize{10}','t_0 = ',num2str(mean(D.refT0{kk}),'%.1f'),' (ns)'])
% title({'\fontsize{14}GPR Snow Properties','\fontsize{12}Density (g/cm^3)',['\fontsize{10}','t_0 = ',num2str(mean(D.refT0{1}),'%.1f'),' (ns)'] })
set(gca,'fontsize',10,'fontweight','bold','xtick',...
    [mean(D.refRho{kk})-1.96.*std(D.refRho{kk}),mean(D.refRho{kk}),mean(D.refRho{kk})+1.96.*std(D.refRho{kk})])
xticklabels(round([mean(D.refRho{kk})-1.96.*std(D.refRho{kk}),mean(D.refRho{kk}),mean(D.refRho{kk})+1.96.*std(D.refRho{kk})],2))
ylim([0 0.35])
grid on
end
end