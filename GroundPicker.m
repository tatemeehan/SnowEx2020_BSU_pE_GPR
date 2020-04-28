%% Automatic Ground Picker
deconIx = 10;
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
%     coPick = tmpCoherence(groundIx);
    % Save Picks
    D.groundIx{ii} = groundPick;
    D.Time2Ground{ii} = groundPick.*D.dt;
end
% Noise Region Picker
% for ii = 1:MD.nFiles
%     topPick = zeros(size(D.Coherence{ii},2));
%     for jj = 1:size(D.Coherence{ii},2)
%         % Adaptive Treshold
%         trc = D.Coherence{ii}(:,jj);
%         threshold = quantile(trc,.99);
%         pickIx = find(trc > threshold);
% %         filtPickIx = find(pickIx > quantile(pickIx,.05) & pickIx < quantile(pickIx,0.95));
%         % Cluster Weighting
% %         if std(pickIx) > 1/D.dt
%         groundPick(jj) = mean(pickIx);
%     end
%     D.
%     D.Time2Ground = groundPick.*D.dt;
% end
        % Cluster Weighting
%         if std(pickIx) > 1/D.dt
%         filtPickIx = find(pickIx > quantile(pickIx,.05) & pickIx < quantile(pickIx,0.95));

 % Filter Picks using Savitzky-Golay Filter
%     ord = 3;
%     R = 25;
%     alpha = 30; % Attenuation dB
%     if alpha > 50
%     beta = 0.1102.*(alpha-8.7);
%     elseif alpha <= 50 && alpha>= 21
%         beta = 0.5842.*(alpha-21).^(0.4)+0.07886.*(alpha-21);
%     else
%         beta = 0;
%     end
%     w = kaiser(R,beta);
%     groundPickSG = sgolayfilt(groundPick,ord,R,w);

    clear('tmpCoherence','groundPick','threshold','pickIx','trc','topPick','deconIx','coPick','R')
