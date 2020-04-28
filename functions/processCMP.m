function [ Rad, TWT ] = processCMP( Rad, f0, dt, TWT, offset )
% processCMP is a data processing Subroutine
%   
%   Written by Tate Meehan, Boise State University, NASA SnowEx 2020
%
% Process Data
% Filter Parameters are Established within Function
isDisplay = 0;       % Control Text Display to Screen
isNormalize = 0;     % Control Flag to Normalize Data
isMedFilt = 0;       % Control Flag for Median Subtraction Filter
isBandPass = 1;      % Control Flag to Band-Pass Filter Data
isTimeZero = 1;      % Control Flag for Time-Zero Correction
isSVDfilt = 0;       % Control Flag for SVD Component Subtraction Filter
isSpatialMedFilt = 0;% Control Flag for Spatial Median Subtraction Filter
isExpGain = 1;       % Control Flag for Ramped Gain of Data
isAGCgain = 0;       % Control Fla for AGC Gain of Data
isStak = 0;          % Control Flag to Stack Data

    %----------------------------------------------------------------------      
    % Convert Units
    f0Hz = f0 * 1e6;        % [Hz]
    dtSec = dt * 1e-9;      % [s]
    [nsamp, ntrcs] = size(Rad);

    %----------------------------------------------------------------------
    % Normalize Data
    if isNormalize
        if isDisplay
        display( 'Begin Normalize')
        tic
        end
        
        Rad = normalize( Rad );
        
        if isDisplay
        display( 'Normalize Done')
        toc
        display(' ')
        end
    end    
    
        
    %----------------------------------------------------------------------
    % Median Subtraction Filter
    if isMedFilt
        if isDisplay
        display( 'Begin Median Subtraction Filter')
        tic
        end
        
        % Parameters
        % Rank of Median Subtraction Filter
        % Nominal Frequency Pass
        MedFiltR = 2.*(ceil(1/((f0Hz)*dtSec))-1)+1;
        % Low Pass
%          MedFiltR = 2.*(ceil(1/((f0Hz/2.5)*dtSec))-1)+1;
        % High Pass
%           MedFiltR = 2.*(ceil(1/((f0Hz*2.5)*dtSec))-1)+1;

        % Function
        Rad = medfilt1( Rad, MedFiltR, [], 2,'omitnan','truncate' );
        
        if isDisplay
        display( 'Median Subtraction Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------    
    % Band-Pass Filter
    if isBandPass
        if isDisplay
        display( 'Begin Band-Pass Filter')
        tic
        end
        
        % Parameters
        % Build Two Octave Filter About Nominal Frequency
        fMin = f0Hz/2; % [Hz]   
        fMax = f0Hz*2; % [Hz]
%         fMin = f0Hz/2; % [Hz]   
%         fMax = f0Hz; % [Hz]
        
        % Function
        Rad = bpfilter( Rad, dtSec, fMin, fMax, 8 );
        
        if isDisplay
        display( 'Band-Pass Filter Done')
        toc
        display(' ')
        end
    end
            
    %----------------------------------------------------------------------
    % Time-Zero Correction
    if isTimeZero
        if isDisplay
        display( 'Begin Time Zero Correction')
        tic
        end
        
        
        % Parameters      
        merR = 1;   %Rank of MER window [ns]
        powMER = 3; % Power of MER operation
        
       
        % Function
        [Rad, TWT, ~] = timeZero( Rad, TWT, dt, merR, powMER, offset, 1);        
        
        if isDisplay
        display( 'Time Zero Correction Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % SVD Principle Component Subtraction Filter
    if isSVDfilt
        if isDisplay
        display( 'Begin Singular Value Decomposition Filter')
        tic
        end
        PCA = 0.75; % PCA percentage threshold
        Rad = SVDSfilter( Rad, PCA );
        
        if isDisplay
        display( 'Singular Value Decomposition Filter Done')
        toc
        display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Spatial Median Subtraction
    %     SpatialMedianSubtract
    if isSpatialMedFilt
        display('Begin Coherent Noise Subtraction')
        tic
        [ Rad ] = movingMedianSubtraction( Rad, round(size(Rad,2).*.75) );
       
        display('Coherent Noise Removal Done')
        toc
        display(' ')
    end


    %----------------------------------------------------------------------
    % Exponential Time Dependant Gain
    if isExpGain
        if isDisplay
        display( 'Begin Power Gain')
        tic
        end
        
        % Parameters
        tpow = 2; % Filter Order, 1 is Linear
        rollOff = 75;
        
        % Function
        Rad = gain(Rad, tpow, rollOff );
        
        if isDisplay
        display( 'Power Gain Done')
        toc
        display(' ')
        end
    end
    
    % ---------------------------------------------------------------------
    % Automatic Gain Control
    if isAGCgain
        if isDisplay
            tic
            display( 'Begin AGC')
        end
        
        % Parameters
        R = ceil(nsamp/2);
        type = 2; % Trace Normalize: 0 = None, 1 = amplitude, 2 = RMSnorm

        Rad = AGCgain(Rad,R,type);
        
        if isDisplay
            display( 'AGC Done')
            toc
            display(' ')
        end
    end
    
    %----------------------------------------------------------------------
    % Stacking Filter
    if isStak
        if isDisplay
        display( 'Begin Stack')
        tic
        end
        
        % Parameters
        StakFiltR = 10; % Filter Rank
        
        % Function
        Rad = StakR(Rad,StakFiltR);
        
        if isDisplay
        display( 'Stacking Done')
        toc
        display(' ')
        end
    end

end

