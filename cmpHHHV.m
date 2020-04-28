%% Common Midpoint Gather Dual-Pole (HH-HV) GPR (cmpHHHV)
clear; close all; clc;
%% Meta Data is a structure MD
MD.dataDir = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR\02012020\CMP1';
% MD.dataDir = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR\02012020\CMP2';
% MD.dataDir = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR\02012020\CMP3';

% Import Utilities
addpath './functions';
addpath './colormaps';
isSave = 0;

% Extract Files
MD.workDir = pwd;
MD.fileNames = dir([MD.dataDir,'/','*.nc']);
MD.lineNo = 1:length(MD.fileNames);    % Array of data "LINE" numbers
MD.nFiles = length(MD.lineNo);        % Number of Files

% Establish Tx Rx Geometry for CMP Gathering
MD.nTx = 1;               % Number of Transmitters in Sequence
MD.nRx = 2;               % Number of Receivers in Sequence

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
    % Choose Polarization
    polAns = questdlg('Choose the CMP Polarization.','CMP Polarization','HH','HV','HH');
    MD.polIx = find(strcmp({polAns},MD.Polarization));
    % Choose Number of Horzons to Pick
    prompt = {'Number of Direct Wave Horizons','Numer of Reflection Horizons','Direct Wave Norm', 'Reflection Norm' };
    dlgtitle = 'Input Horizons';
    dfaultinput = {'1','1','L2','L2'};
    dim = 1;
    horizons = inputdlg(prompt,dlgtitle,dim,dfaultinput);
    MD.ndir = str2num(horizons{1}); % Number of Direct Horizons
    MD.nref = str2num(horizons{2}); % Number of Reflection Horizons
    % Choose Inversion norm
    MD.dirNorm = horizons{3};
    MD.refNorm = horizons{4};
    
    try close(998)
        close(998)
    catch
    end
    % Pick Direct Wave
    for kk = 1:MD.ndir
        picks = ManualPicker(D.CMP{MD.polIx},'zero',3,MD.dt,MD.offsetArray,MD.TimeAxis);
        D.dirTpick{kk} = picks;
    end
    
    % Pick Ground Reflection
    for kk = 1:MD.nref
        picks = ManualPicker(D.CMP{MD.polIx},'zero',3,MD.dt,MD.offsetArray,MD.TimeAxis);
        D.refTpick{kk} = picks;
    end
    
    clear('dim','dfaultinput','prompt','polAns','horizons','picks')
    pause(5)
    close(998)
    %% BootStrapping Velocity Analysis
    MCix = 1:length(MD.offsetArray);
    nMC = 250; % Number of MC Simulations
    for kk = 1:nMC
        binIx = datasample(MCix,round(0.5.*length(MD.offsetArray)));
        % Least Squares Velocity Analysis
        % Direct Wave Analysis
        for hh = 1:MD.ndir
            G = [ones(length(MD.offsetArray(binIx)),1),MD.offsetArray(binIx)'];
            d = D.dirTpick{hh}(binIx);
            if strcmp(MD.dirNorm,'L2')
                % L2 norm
                m = G\d;
            else
                %L1 norm
                p = 1; % p-Norm
                tolr = 1e-2; % Residual Tolerance
                tolx = 5e-4; % Convergence Tolerance
                maxiter = 100; % Maximm Iterations
                m = irls(G,d,tolr,tolx,p,maxiter);
            end
            % Save Intercept Time, Velocity, Wavelength, and Density
            D.dirT0{hh}(kk) = m(1);
            D.dirV{hh}(kk) = 1/m(2);
            D.dirZ{hh}(kk) = D.dirV{hh}(kk)./MD.f0GHz;
            D.dirRho{hh}(kk) = DryCrim(D.dirV{hh}(kk));
        end
        % Reflecttion Velocity Analysis
        for hh = 1:MD.nref
            G = [ones(length(binIx),1),(MD.offsetArray(binIx))'.^2];
            d = [(D.refTpick{hh}(binIx) - D.dirT0{1}(kk)).^2];
            if strcmp(MD.refNorm,'L2')
                % L2 norm
                m = G\d;
            else
                % L1 norm
                p = 1; % p-Norm
                tolr = 1e-2; % Residual Tolerance
                tolx = 5e-4; % Convergence Tolerance
                maxiter = 100; % Maximm Iterations
                m = irls(G,d,tolr,tolx,p,maxiter);
            end
            % Save Intercept Time, Velocity, Depth, and Density
            D.refT0{hh}(kk) = sqrt(m(1));
            D.refV{hh}(kk) = 1./sqrt(m(2));
            D.refZ{hh}(kk) = D.refT0{hh}(kk).*D.refV{hh}(kk)./2;
            D.refRho{hh}(kk) = DryCrim(D.refV{hh}(kk));
        end
    end
    clear('G','d','m','ix','picks','binIx')
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