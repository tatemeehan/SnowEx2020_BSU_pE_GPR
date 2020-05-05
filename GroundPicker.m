%% Automatic Ground Picker
deconIx = 10; % Subtract integer samples from pick (layman's deconvolution)
for ii = 1:MD.nFiles
    groundPick = zeros(size(D.Coherence{ii},2),1);
    tmpCoherence = D.Coherence{ii};
    for jj = 1:size(D.Coherence{ii},2)
        % Adaptive Treshold
        trc = tmpCoherence(:,jj);
        threshold = quantile(trc,.99);
        pickIx = find(trc > threshold);
        groundPick(jj) = (mean(pickIx)-deconIx);
    end
    % Median Filter 
    R = 25;
    groundPick = ceil(medfilt1(groundPick,R));
    groundIx = sub2ind(size(tmpCoherence),groundPick,[1:size(tmpCoherence,2)]');
    % Save Picks
    D.groundIx{ii} = groundPick;
    D.Time2Ground{ii} = groundPick.*D.dt;
end
    clear('tmpCoherence','groundPick','threshold','pickIx','trc','topPick','deconIx','coPick','R')
