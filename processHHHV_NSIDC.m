%% Read Sensors and Software Data
% This Script Reads the Binary Multi-channel GPR Data and Performs the
% pre-processing and signal processing routines. 
% GPS Locations are incorporated in NSIDC Sea Ice Polarstereographic North
%
% The Array Geometry is Installed
% The Acquisition Method is Decided
% The Near Channel (7) is Killed
% The Time Window can be trimmed
% Data is de-multiplex (grouped into channels) and coarse trace shifts
% are applied as the first cut at correcting for digital time-sampling 
% errors.
% Inside the function wrapper processCommonOffset.m de-WOW filtering and
% trace stacking are applied. There are Various filter parameters decided
% in this function.

    D.Rad = cell(1,MD.nFiles);
    D.trhd = cell(1,MD.nFiles);
%     GPS = cell(1,MD.nFiles);
%     Year = cell(1,MD.nFiles);
    TimeAxis = cell(1,MD.nFiles);
    
    for ii = 1 : MD.nFiles
        tic
        %------------------------------------------------------------------
        % Multiplexed Channel Record
        filename = MD.fileNames(ii).name;
        filepath = fullfile(MD.dataDir,filename);
        % Read netCDF data file
        disp(' ')
        fprintf('Reading .nc File \n')
        tic
        ncRad = ncread(filepath,'DATA');
        D.trhd{ii} = ncRad(1:20,:);
        D.Rad{ii} = ncRad(21:end,:);
        clear('ncRad');
        D.f0 = (D.trhd{ii}(3,1)); % [MHz]
        D.dt = D.trhd{ii}(4,1); % [ns]
%         D.dx = diff(D.trhd{ii}(2,1:2)); % [m]
        % Need to Automate Offset Array from .nc File
%         offsetArray = unique(D.trhd{ii}(22,1:100));

 
        % Nominal Frequency GHz
        D.f0GHz = D.f0/1000;
        % No. Traces in Multiplexed Data
        [~, multiplexNtrcs] = size(D.Rad{ii});

        % Input Polaraztion Configureation
        % SkiPol
        D.Polarization = {'HH','HV'};

        
        D.offsetArray = D.trhd{ii}(7,chan);      % Offset Array for CMP Gather
        nChan = length(D.offsetArray);           % Refresh nChan before kill
        
        % Pad Data with Instrument Zero
        padding = 0;
        instrumentPad = zeros(padding,size(D.Rad{ii},2));
        if padding ~= 0
            for jj = 1:size(D.Rad{ii},2)
                instrumentZero = D.Rad{ii}(1,jj);
                instrumentPad(:,jj) = ones(padding,1).*instrumentZero;
            end
            D.Rad{ii} = [instrumentPad;D.Rad{ii}];
        end  
        
        % Trim Time Window
        if isTrimTWT
            endSamp = 200; % The new end of the columns
            reSample = endSamp + padding;   % Number of Wanted Samples
            D.Rad{ii} = D.Rad{ii}(1:reSample,:);
        end
        
        % Allocation Here
        if ii == 1
            Radar = cell(nChan,MD.nFiles); traceIx = cell(nChan,MD.nFiles);
%             Array = cell(nChan,MD.nFiles);
        end
        
        for jj = chan
            % DeMux Sequential Data
            % GPS DeadReckoning completed in preProcessing
            [Radar{jj,ii},D.trhd{ii},traceIx{jj,ii},D.Distance{jj,ii}] = deMuxNSIDC(D.Rad{ii},D.trhd{ii},chan(jj));
        end
        
        % Remove Every Nth Trace for Data Reduction
        rmNtrc = 2;
        dt = D.dt;
%         parfor (jj =  1:nChan, nWorkers)
        for jj = 1:nChan
                % Reduce Data Volume
                if isReduceData
                    Radar{jj,ii} = Radar{jj,ii}(:,1:rmNtrc:end);
                    traceIx{jj,ii} = traceIx{jj,ii}(1:rmNtrc:end);
                end
                
                % Process Common Offset Channels
%                 disp(' ')           
%                 fprintf(['Begin Signal Processing in Common-Offset Domain ',...
%                     filename, ' CHAN0', num2str(jj),'\n'])
            % Store Travel-Time Axis
            TWT = [0:dt:(dt.*(size(Radar{jj,ii},1)-1))]';

            [Radar{jj,ii}, TimeAxis{jj,ii}] = processCommonOffset(Radar{jj,ii}, ...
                D.f0, D.dt, TWT, D.offsetArray(jj) );

        end
        if isReduceData
            tmpIx = sort(cat(2,traceIx{:,ii}));
            D.trhd{ii} = D.trhd{ii}(:,tmpIx);
        end
        
        % Trim all time windows to same length
        for jj = chan
            if jj == 1
            minIx = length(TimeAxis{jj,ii});
            minChan = jj;
            elseif minIx > length(TimeAxis{jj,ii})
                minIx = length(TimeAxis{jj,ii});
                minChan = jj;
            end
        end
        clear('TimeAxis')
        % Store Travel-Time Axis
            TimeAxis = [0:dt:(dt.*(minIx-1))]';
        % Trim Channels to Consistent Travel Time Axis
        for jj = chan
            Radar{jj,ii} = Radar{jj,ii}(1:minIx,:);            
        end
        D.Radar = Radar;
        kk = 0;
        D.X{ii} = zeros(length(1:nChan:size(D.trhd{ii},2)),1);
        D.Y{ii} = D.X{ii}; D.Z{ii} = D.X{ii};D.DistanceAxis{ii}=D.X{ii};
        % Extract Positions
        for jj = 1:nChan:size(D.trhd{ii},2)
            kk = kk + 1;
                D.X{ii}(kk) = D.trhd{ii}(9,jj);
                D.Y{ii}(kk) = D.trhd{ii}(10,jj);
                D.Z{ii}(kk) = D.trhd{ii}(11,jj);
                D.DistanceAxis{ii}(kk) = D.trhd{ii}(12,jj);
        end
        D.TimeAxis{ii} = TimeAxis;
        clear('TimeAxis')
        
        fprintf('Signal Processing Done \n')
        toc
        display(' ')
    end
    MD.nChan = nChan; MD.chan = chan;
    D = rmfield(D,'Distance');
    clear('Rad','tmpIx','Distance','endDist','chan','nChan','TWT','dt','TimeAxis','minIx','minChan',...
        'dH','dS','dT','dupIx','dX','dY','dZ','gatherLength','instrumentPad',...
        'padding','midPointArray','multiplexNtrcs','pol','polOffset','Radar','rmNtrc',...
        'rxGeo','spaceIx','spaceTime','traceIx','traceMod','txGeo','X','xArray','xTrc','Y','Z','vThreshold');