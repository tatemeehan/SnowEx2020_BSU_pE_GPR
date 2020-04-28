function [dat,hdr1,trhd,dt,f0,xArray,dx,date] = readSensorsSoftwareData(filename)
%
% USAGE: [dat,hdr1,trhd,dt,f0,xArray] = readSensorsSoftwareData(filename)
%
% Description: Read a sensors and software GPR file, both hdr and dt1
% components. The data matrix is returned as well as the header
% information.
%
% INPUT:
%   filename = full filename including the path if not in current directory
% OUTPUT:
%   dat    = data matrix ( row=time, column=space)
%   hdr1   = survey geometry headers
%   trhd   = trace headers
%   dt     = sample interval [ns]
%   f0     = antenna central frequency
%   xArray = position vector corresponding to columns in dat
%
% Written by: Dylan Mikesell (dylanmikesell@boisestate.edu)
% Last modified: 2 November 2016
% By: Tate Meehan - dt calculation was edited
% Original: dt  = t / ( ns - 1 );
% Modified: dt  = t / ( ns );
% Based on code from John Bradford
%
% V1: wrote code to read 500 MHz demultiplexed data
% V2: rewrote code to handle non-demultiplexed 100 MHz data as well. This
% is a more robust function now but assumes that headers are the same
% within the two data types.

%--------------------------------------------------------------------------
% Load the HEADER file
hdfile = [filename '.HD'];

if ~exist(hdfile,'file')
    disp(hdfile)
    error('MATLAB/readSensorsSoftwareData: HEADER file does not exist.');
end

fid = fopen(hdfile,'r');
hdr1 = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
hdr1 = hdr1{1}; % a cell array of strings, one per line in the file.

% Find all the blank lines using cellfun and remove
hdr1( cellfun( @isempty, hdr1 ) ) = [];
% Grab Date Stamp
date = hdr1{3};
% remove the first 3 header lines
hdr1(1:3) = [];

for ii = 1 : numel(hdr1) 
    C{ii} = strtrim( strsplit(hdr1{ii},'=') ); % separate variables and trim whitespace
end

% Find the parts we need
ntr = round( str2double( C{1}{1,2} ) ); % number of traces 
ns  = round( str2double( C{2}{1,2} ) ); % number of sample per trace 
t   = str2double( C{4}{1,2} ); % [nanosecond] t_max 
dt  = t / ( ns ); % [nanosecond] time sample interval

xstart = str2double( C{5}{1,2} ); % [m] start of survey
xstop  = str2double( C{6}{1,2} ); % [m] end of survey
dx     = str2double( C{7}{1,2} ); % [m] sample interval in space
xArray = xstart : dx : xstop; % [m] position vector

f0 = round( str2double( C{9}{1,2} ) ); % [MHz] frequency 

%--------------------------------------------------------------------------
% Load the DATA file
dtfile = [filename '.DT1'];

if ~exist(dtfile,'file')
    disp(dtfile)
    error('MATLAB/readSensorsSoftwareData: DATA file does not exist.');
end

% allocate
dat  = zeros( ns, ntr );
trhd = zeros( 25, ntr );

fid=fopen(dtfile,'r','n');
k = 0;
while (k < ntr)
    trhd(1:22,k+1) = fread( fid, 22, 'real*4' ); % read trace header
    fread( fid, 1, 'integer*1' ); % read blank 
    trhd(23,k+1)    =fread(fid,1,'integer*1'); % read trace header
    fread(fid,2,'integer*1'); % read blank
    trhd(24:25,k+1) = fread(fid,2,'real*4'); % read trace header
    comment         = fread(fid,28); % read trace 'comment'
    dat(:,k+1)      = fread(fid,ns,'int16'); % read data
    k = k + 1; % upate counter
end
fclose(fid);

fprintf('Done reading %s.\n',filename);

return



% Old V1 code
%
% % read the first part of the header
% fid = fopen(hdfile,'r');
% hdr1 = textscan( fid, '%21c%f%*[^\n]', 7, 'headerlines', 3 );
% fclose(fid);
% 
% % extract header values that we want
% S = cell2mat( hdr1(2) );
% 
% ntr = round( S(1) ); % number of traces 
% ns  = round( S(2) ); % number of sample per trace
% t   = S(4); % [nanosecond] t_max
% dt  = t / ( ns - 1 ); % [nanosecond] time sample interval
% 
% % xstart = S(5);
% % xstop = S(6);
% % dx = S(7);
% xArray = S(5) : S(7) : S(6); % [m] position vector
% 
% % read the second part of the header
% fid = fopen(hdfile,'r');
% hdr2 = textscan( fid, '%21c%f%*[^\n]', 4, 'headerlines', 11 );
% fclose(fid);
% 
% % extract header values that we want
% S = cell2mat( hdr2(2) );
% f0 = S( 1 );
