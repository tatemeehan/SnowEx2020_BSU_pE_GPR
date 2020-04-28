%% Multi-channel Dual-Pole (HH-HV) GPR (MxHHHV)
clear; close all; clc;
%% Meta Data is a structure MD
% MD.dataDir = 'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\LDP\122019\GPR_DATA';
% MD.dataDir = 'D:\GrandMesaGPR\nc';
directories = {'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR\01282020',...
    'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR\01292020',...
    'D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR\02042020'};
addpath './functions';
addpath './colormaps';
isSavePicks = 0;
isLoadPicks = 0;
%% Load and concatenate .csv Picks
if isLoadPicks
f1 = readtable([directories{1},'\BSU_pE_GPR_01282020_TWT.csv']);
f2 = readtable([directories{2},'\BSU_pE_GPR_01292020_TWT.csv']);
f3 = readtable([directories{3},'\BSU_pE_GPR_02042020_TWT.csv']);
% Concatenate Tables
f4 = [f1;f2;f3];
% Decimate Data
f5 = f4(1:100:end,:);
% Write the Concatenated File
if isSavePicks
    cd('D:\CRREL_SnowCompaction\CRREL\SnowEx2020\GrandMesa\GPR')
    writetable(f4,'BSU_pE_GPR_01282020_01292020_02042020_TWT.csv')
    writetable(f5,'BSU_pE_GPR_01282020_01292020_02042020_TWT_decimated.csv')
end
end
%% Process Traveltime Data
for hh = 1:size(directories,2)
MD.dataDir = directories{hh};
MD.workDir = pwd;
MD.fileNames = dir([MD.dataDir,'/','*.nc']);
MD.lineNo = 1:length(MD.fileNames);                   % Array of data "LINE" numbers
MD.nFiles = length(MD.lineNo);        % Number of Files
nChan = 2;                      % Number of Recorded Channels
chan =  1:nChan;                % Linear Array of Record Channels

% Establish Tx Rx Geometry for CMP Gathering
MD.nTx = 1;               % Number of Transmitters in Sequence
MD.nRx = 2;               % Number of Receivers in Sequence
%% Processing WorkFlow Controls
% Parallel Computing Enabled
isParallel = 1;
% Read Data
isReadNC = 1;                  % Read Multiplexed Data
isTrimTWT = 1;          % Truncate Recorded Data
isReduceData = 0;

% Write SWE Data
isWrite = 0;

% Load Color Maps
yetBlack = load('yetBlack.txt');
%% Control Parallelization
if isParallel
    % Wake Parallel Computing
    if isempty(gcp('nocreate'))     
        nWorkers = nChan;
        p = parpool(nWorkers);
    else nWorkers = nChan;
    end
else
    nWorkers = 1;
end

%% Read & Process GPR Data
if isReadNC
    processHHHV_NSIDC    
end

%% Measure Ground Reflection Coherency
fprintf('Calculating HH-HV Coherence \n')
tic

QPcoherence

toc
%% Pick Ground Return
fprintf('Picking the Ground \n')
GroundPicker
% Save Picks to .csv
if isSavePicks
    disp('Writing Travel-Time Data')
    tic
% Create Date Frame Header
header = {'Easting','Northing','Elevation','TWT'};
header = strjoin(header,',');
for ii = 1:MD.nFiles
    if ii == 1
    F = [D.X{ii},D.Y{ii},D.Z{ii},D.Time2Ground{ii}];
    else
        % Concatenate to Daily Data Frame
        F = [F;[D.X{ii},D.Y{ii},D.Z{ii},D.Time2Ground{ii}]];
    end
    out = [D.X{ii},D.Y{ii},D.Z{ii},D.Time2Ground{ii}];
    fname = [MD.fileNames(ii).name(1:end-3),'_TWT.csv'];
    cd(MD.dataDir) % Change to Save Directory
    fid = fopen(fname,'w'); % open Filename
    fprintf(fid,'%s\n',header); % Write Headers
    fclose(fid);     % Close File
    % Write Data
    dlmwrite(fname,out,'-append','delimiter',',','precision','%.3f');
    cd(MD.workDir) % Change to Working Directory
    if ii == MD.nFiles
    % Open File and Write Header
    cd(MD.dataDir)
    fname = [MD.fileNames(ii).name(1:end-6),'_TWT.csv'];
    fid = fopen(fname,'w');
    fprintf(fid,'%s\n',header);
    fclose(fid);
    % Write Data
    dlmwrite(fname,F,'-append','delimiter',',','precision','%.3f');
    cd(MD.workDir)
    end
    clear out
end
clear('F','header','date','fname','name','D','MD')
toc
disp(' ')
end
end
%% AGC Gain
isAGCgain = 0;
if isAGCgain
D.RadarAGC = cell(MD.nChan,MD.nFiles);
for ii = 1:MD.nFiles
    for jj = 1:MD.nChan
        D.RadarAGC{jj,ii} = AGCgain(D.Radar{jj,ii},250,2);
    end
end
end
%% Calculate SWE
% pitDataDir = '/home/tatemeehan/git-repository/GrandMesaMarch2019/';
pitDataDir = 'D:\git-repository\GrandMesaMarch2019\';
pitFilename = 'PitSummary.csv';
% Distribute Density Data
pitData = readtable([pitDataDir,pitFilename]);
pitE = pitData.Easting;
pitN = pitData.Northing;
pitDensity = pitData.MeanDensity_kgpm3_;
for ii =  1:MD.nFiles
L = 11;
r = 1000; % 1000 m search radius 
% Averge Density from Snow Pits using Inverse Distance Weighting
D.snowDensity {ii} = distributeDensity(pitE,pitN,pitDensity,D.X{ii},D.Y{ii},r);
D.velocity{ii} = DryCrimVRMS(D.snowDensity{ii});
D.snowDepth{ii} = 100.*((D.Time2Ground{ii}.*D.velocity{ii})./2);
D.SWE{ii} = (D.snowDensity{ii}.*D.snowDepth{ii})./100;
% Smooth SWE
% D.SWE{ii} = movmean(D.snowDensity.*D.snowHeight,L);
end

%% Write DataFrame.csv
if isWrite
    disp('Writing Data Frame')
    tic
% Create Date Frame Header
header = {'Easting','Northing','TWT','Depth','MeanDensity','SWE'};
header = strjoin(header,',');
for ii = 1:MD.nFiles
    if ii == 1
    F = [D.X{ii},D.Y{ii},D.Time2Ground{ii},D.snowDepth{ii},...
        D.snowDensity{ii},D.SWE{ii}];
    else
        % Concatenate to Daily Data Frame
        F = [F;[D.X{ii},D.Y{ii},D.Time2Ground{ii},D.snowDepth{ii},...
            D.snowDensity{ii},D.SWE{ii}]];
    end
    if ii == MD.nFiles
    % Open File and Write Header
    cd(MD.dataDir)
    name = 'PulseEKKO_QP_';%
    date = MD.fileNames(MD.lineNo(ii)+1).name;
    fname = [name,date(11:21),'.csv'];
    fid = fopen(fname,'w');
    fprintf(fid,'%s\n',header);
    fclose(fid);
    % Write Data
    dlmwrite(fname,F,'-append','delimiter',',','precision','%.1f');
    cd(MD.workDir)
    end
end
clear('F','header','date','fname','name')
toc
disp(' ')
end
%% Sum Radar Energy above Ground
isEnergyRatio = 0;
if isEnergyRatio
% Sum energy of Coherency Grams
for ii =  1:MD.nFiles
    groundIx = D.groundIx{ii};
    Co = D.HHHVCoherence{ii};
    Cx = D.VVVHCoherence{ii};
    tmp = zeros(size(D.Radar{1,ii},2),1);
    parfor (kk = 1:size(D.Radar{1,ii},2),nWorkers)
        tmp(kk) = sum(Cx(1:(groundIx(kk)-deconIx),kk))./sum(Co(1:(groundIx(kk)-deconIx),kk));
    end
    D.EnergyRatio{ii} = movmean(tmp,5.*L);
end
tic
energySum = cell(MD.nChan,MD.nFiles);
EnergyRatio = cell(MD.nChan,MD.nFiles);
Rad= D.Radar;
rmIx = 25; % Remove 25 samples above coherent reflector
for ii = 1:MD.nFiles
%     rad = D.Radar{jj,ii};
    for jj = 1:MD.nChan
        for kk = 1:size(D.Radar{jj,ii},2)
%             energySum{jj,ii}(kk) = sum((Rad{jj,ii}([1:D.groundIx(:)-rmIx],kk).^2));
            energySum{jj,ii}(kk) = sum((Rad{jj,ii}(:,kk).^2));
        end
    end
    for jj = 1:MD.nChan
        if strcmp(D.Polarization{ii,jj},'HH')
            % Cross-Pole HV to Co-Pole HH
            EnergyRatio{jj,ii} = energySum{1,ii}./energySum{jj,ii};
        elseif strcmp(D.Polarization{ii,jj},'HV')
            % HV to VV
            EnergyRatio{jj,ii} = energySum{jj,ii}./energySum{3,ii};
        elseif strcmp(D.Polarization{ii,jj},'VV')
            % VH to VV
            EnergyRatio{jj,ii} = energySum{4,ii}./energySum{jj,ii};
        elseif strcmp(D.Polarization{ii,jj},'VH')
            % VH to HH
            EnergyRatio{jj,ii} = energySum{jj,ii}./energySum{2,ii};
        else
            error('Polarizations are not Defined!')
        end 
    end
end
toc
end
%% Make Figures
isMakeFigures = 0;
if isMakeFigures
    set(0,'DefaultAxesFontName','Serif')
set(0,'DefaultTextFontName','Serif')
% set(0,'DefaultAxesFontName','FreeSerif')
% set(0,'DefaultTextFontName','FreeSerif')
% set(0,'DefaultFontSize',14,'DefaultFontWeight','bold')
xix = [find(D.DistanceAxis{1}>3700,1):find(D.DistanceAxis{1}>5300,1)];
% yix = [1:300;
% xix = 1:length(D.DistanceAxis{1});
figure();
pcolor(D.DistanceAxis{ii}(xix),D.TimeAxis{ii},D.Coherence{ii}(:,xix));
shading interp; colormap(yetBlack);axis ij;hold on;caxis([0,1]);
% imagesc(D.DistanceAxis(xix),D.TimeAxis,D.Coherence{1}(:,xix));
% shading interp; colormap(yetBlack);hold on;
plot(D.DistanceAxis{ii}(xix),D.Time2Ground{ii}(xix),'.','color',[.7,.7,.7])
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Quad-Pol Coherence';
c.FontSize = 12; c.Label.FontSize = 14;
set(gca,'fontsize',14,'fontweight','bold','layer','top')
%          title('Coherence')
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
set(gca,'fontsize',14,'fontweight','bold')

figure();
pcolor(D.DistanceAxis(xix),D.TimeAxis,D.HHHVCoherence{1}(:,xix));
shading interp; colormap(yetBlack);axis ij;hold on;caxis([0,1]);
plot(D.DistanceAxis(xix),D.Time2Ground(xix)-1,'.','color',[.7,.7,.7])
c = colorbar; c.Location = 'northoutside';c.Label.String = 'HH-HV Coherence';
c.FontSize = 12; c.Label.FontSize = 14;
set(gca,'fontsize',14,'fontweight','bold','layer','top')
%          title('Coherence')
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
set(gca,'fontsize',14,'fontweight','bold')

figure();
pcolor(D.DistanceAxis(xix),D.TimeAxis,D.VVVHCoherence{1}(:,xix));
shading interp; colormap(yetBlack);axis ij;hold on;
plot(D.DistanceAxis(xix),D.Time2Ground(xix)-2.5,'.','color',[.7,.7,.7])
c = colorbar; c.Location = 'northoutside';c.Label.String = 'Cross-Pol Coherence';
c.FontSize = 12; c.Label.FontSize = 14;
set(gca,'fontsize',14,'fontweight','bold','layer','top')
%          title('Coherence')
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
set(gca,'fontsize',14,'fontweight','bold')

% Plot Time Gain RadarGram HV
cq = quantile(D.Radar{1}(:),[0.005,0.995]);
figure();
pcolor(D.DistanceAxis{ii}(xix),D.TimeAxis{ii},D.Radar{1,1}(:,xix));
shading interp;colormap(bone);axis ij;
caxis(cq);
set(gca,'fontsize',14,'fontweight','bold','layer','top')
title('HV Polarization')
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
set(gca,'fontsize',14,'fontweight','bold')

% Plot Time Gain RadarGram HH
cq = quantile(D.Radar{2}(:),[0.005,0.995]);
figure();
pcolor(D.DistanceAxis{ii}(xix),D.TimeAxis{ii},D.Radar{2,1}(:,xix));
shading interp;colormap(bone);axis ij;
caxis(cq);
set(gca,'fontsize',14,'fontweight','bold','layer','top')
title('HH Polarization')
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
set(gca,'fontsize',14,'fontweight','bold')

% Plot Time Gain RadarGram VV
cq = quantile(D.Radar{3}(:),[0.005,0.995]);
figure();
pcolor(D.DistanceAxis,D.TimeAxis,D.Radar{3,1});
shading interp;colormap(bone);axis ij;
caxis(cq);
set(gca,'fontsize',14,'fontweight','bold','layer','top')
title('VV Polarization')
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
set(gca,'fontsize',14,'fontweight','bold')

% Plot Time Gain RadarGram VH
cq = quantile(D.Radar{4}(:),[0.005,0.995]);
figure();
pcolor(D.DistanceAxis,D.TimeAxis,D.Radar{4,1});
shading interp;colormap(bone);axis ij;
caxis(cq);
set(gca,'fontsize',14,'fontweight','bold','layer','top')
title('VH Polarization')
ylabel('Travel-time (ns)')
xlabel('Distance (m)')
set(gca,'fontsize',14,'fontweight','bold')

% Thin Results
corrIx = find(D.EnergyRatio{1}(xix)>quantile(D.EnergyRatio{1}(xix),.05) & D.EnergyRatio{1}(xix)<quantile(D.EnergyRatio{1}(xix),.95));
figure();
scatter(D.EnergyRatio{1}(xix(corrIx)),D.SWE(xix(corrIx)),'k'); hold on;
lsline;
[R,P] = corrcoef(D.EnergyRatio{1}(xix(corrIx)),D.SWE(xix(corrIx)));
annotation('textbox',[.8,.8,.1,.1],'string',{['R^2: ',num2str(R(1,2).^2)]},'FitBoxToText','on');
% ;['p: ',num2str(P)]

figure();plot(D.DistanceAxis(xix),D.EnergyRatio{1}(xix),'.k')
hold on; plot(D.DistanceAxis(xix),D.SWE(xix),'.r')
set(gca,'fontsize',14,'fontweight','bold')
title('HV/HH Energy Ratio and SWE')
xlabel('Distance (m)')
legend('HV/HH','SWE','location','NorthWest')

figure();
scatter(EnergyRatio{1,1}(xix),D.SWE(xix),'k');hold on; lsline
[R,P] = corrcoef(EnergyRatio{1,1}(xix),D.SWE(xix));
annotation('textbox',[.7,.2,.1,.1],'string',{['R^2: ',num2str(R(1,2).^2)]},'FitBoxToText','on','LineStyle','none','fontsize',14);
title('HV/VV')
ylabel('SWE')
xlabel('Energy Ratio (HV/VV)')
set(gca,'fontsize',14,'fontweight','bold')

figure();
scatter(EnergyRatio{2,1}(xix),D.SWE(xix),'k'); hold on; lsline
[R,P] = corrcoef(EnergyRatio{2,1}(xix),D.SWE(xix));
annotation('textbox',[.7,.2,.1,.1],'string',{['R^2: ',num2str(R(1,2).^2)]},'FitBoxToText','on','LineStyle','none','fontsize',14);
title('HV/HH')
ylabel('SWE')
xlabel('Energy Ratio (HV/HH)')
set(gca,'fontsize',14,'fontweight','bold')


figure();
scatter(EnergyRatio{3,1}(xix),D.SWE(xix),'k'); hold on; lsline
[R,P] = corrcoef(EnergyRatio{3,1}(xix),D.SWE(xix));
annotation('textbox',[.7,.2,.1,.1],'string',{['R^2: ',num2str(R(1,2).^2)]},'FitBoxToText','on','LineStyle','none','fontsize',14);

title('VH/VV')
ylabel('SWE')
xlabel('Energy Ratio (VH/VV)')
set(gca,'fontsize',14,'fontweight','bold')

figure();
scatter(EnergyRatio{4,1}(xix),D.SWE(xix),'k');hold on;lsline
[R,P] = corrcoef(EnergyRatio{4,1}(xix),D.SWE(xix));
annotation('textbox',[.7,.2,.1,.1],'string',{['R^2: ',num2str(R(1,2).^2)]},'FitBoxToText','on','LineStyle','none','fontsize',14);
title('VH/HH')
ylabel('SWE')
xlabel('Energy Ratio (VH/HH)')
set(gca,'fontsize',14,'fontweight','bold')

end