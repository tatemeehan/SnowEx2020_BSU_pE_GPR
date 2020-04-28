function [ Rad, trhd, trcIx, xArray ] = deMuxNSIDC( Rad, trhd, chan)
% DeMux - Demultiplexes Sensors and Software Multi-offset Radar Data.
% DeadReckon - Calculates Truer GPR Antenna Positions by High Accuracy GPS. 
  % Assumtions: The GPS Antenna is located at Tx1
  %             The GPR Array is towed straight
  %             The Shot Gather Bin Center is the Array mean
  % Errors caused by these assumtions are significant, however, inevitable.
  % The dead reckon approach is utilized to distribute the GPS position to
  % each GPR antenna.
  
% Boise State University: Tate Meehan, NASA ISGC 2019

    % DeMux Data
        trcIx = find(trhd(2,:)==chan);
        Rad=Rad(:,trcIx);
        xArray = trhd(2,trcIx);
end