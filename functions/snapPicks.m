function [tPick] = snapPicks(signal,rawPick,ixRange,dt,window,method)
%
%       signal = just like it says
%       ixRange  = [samples] range -+ pick
%              dt = [ns] sample rate
%          method = pick the min or max or zero amplitude in window

% Snap Picks to maximum, minimum, or zero amplitude
switch method
    case 'max'
        [~,ix]  = max(signal(floor(ixRange)));
        tPick  = floor(ixRange(ix));
        
    case 'min'
        [~,ix]  = min(signal(floor(ixRange)));
        tPick  = floor(ixRange(ix));
        
    case 'zero'    
        tmpIx1 = find(diff(sign(signal(floor(ixRange)))))+1;
        % Expanding Window Search
        while isempty(tmpIx1)
            window = window + 1;
            ixRange = int8((round(rawPick,1)-window.*dt)./dt):...
                int8((round(rawPick,1)+window.*dt)./dt);
            tmpIx1 = find(diff(sign(signal(floor(ixRange)))))+1;
        end
        % Multiple Zero Crossing in Window Grab Zero Nearest to Cursor Pick
        [~,tmpIx3] = min(abs(floor(double(ixRange(tmpIx1))).*dt - rawPick));
        tmpIx1 = tmpIx1(tmpIx3);
        % Grab Sample Nearest to zero
        [~, tmpIx2] = min(abs(signal(ixRange([tmpIx1-1,tmpIx1]))));
        if tmpIx2 == 1
            ix = tmpIx1-1;
        else
            ix = tmpIx1;
        end
        tPick  = floor(ixRange(ix));   
        
end
% Convert Indicies to Travel-time
tPick = double(tPick).*dt;
end
        
