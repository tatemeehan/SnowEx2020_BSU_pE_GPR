%% Multi-channel Dual-Pole (HH-HV) GPR (MxHHHV)
clear; close all; clc;
%% Meta Data is a structure MD
directories = {'Enter the Directory Containing the Data',...
    'You may enter multiple directories','like this'};
addpath './functions';
addpath './colormaps';
isSavePicks = 0;

%% Process Traveltime Data
for hh = 1:size(directories,2)
MD.dataDir = directories{hh};
MD.workDir = pwd;
MD.fileNames = dir([MD.dataDir,'/','*.nc']);
MD.lineNo = 1:length(MD.fileNames);                   % Array of data "LINE" numbers
MD.nFiles = 1;%;length(MD.lineNo);        % Number of Files
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
isTrimTWT = 1;                 % Truncate Recorded Data
isReduceData = 0;              % Thin Traces

% Write Travel-time Data
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