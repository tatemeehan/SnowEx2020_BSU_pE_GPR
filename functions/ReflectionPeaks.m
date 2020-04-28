function [twt,GrdPick,SurfPick] = ReflectionPeaks(rawsignal,GrdRange,SurfRange,dt,method)
%
%       rawsignal = just like it says
%       GrdRange  = [samples] range -+ pick
%              dt = [ns] sample rate
%          method = pick the min or max or zero amplitude in window

switch method
    case 'max'
        [~,IxG]  = max(rawsignal(floor(GrdRange)));
        GrdPick  = floor(GrdRange(IxG));
        
        [~,IxS]  = max(rawsignal(floor(SurfRange)));
        SurfPick = floor(SurfRange(IxS));
        
    case 'min'
        [~,IxG]  = min(rawsignal(floor(GrdRange)));
        GrdPick  = floor(GrdRange(IxG));
        
        [~,IxS]  = min(rawsignal(floor(SurfRange)));
        SurfPick = floor(SurfRange(IxS));
        
    case 'zero'    
        [IxG]  = knnsearch(rawsignal(floor(GrdRange)),0);
        GrdPick  = floor(GrdRange(IxG));
        
        [IxS]  = knnsearch(rawsignal(floor(SurfRange)),0);
        SurfPick = floor(SurfRange(IxS));
        
end

twt = (GrdPick-SurfPick)*dt;
        
