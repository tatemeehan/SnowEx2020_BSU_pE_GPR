%% Ground Coherence
L = 11; % Window Length
H = hamming(L).*ones(L,MD.nChan);
D.Energy = cell(MD.nFiles,MD.nChan);
isCrossCoherence = 0;

for ii = 1:MD.nFiles
    C = zeros(size(D.Radar{1,ii}));
    Co = C; Cx = C;
    for (kk = 1:size(D.Radar{1,ii},2))
        iter = 1;
        coIter = 1;
        cxIter = 1;
        for jj = 1:MD.nChan
            % Extract Shot Gather
            s(:,iter) = D.Radar{jj,ii}(:,kk);  
            if isCrossCoherence
            if strcmp(D.Polarization{jj},'HH') || strcmp(D.Polarization{jj},'HV')
                sCo(:,coIter) = D.Radar{jj,ii}(:,kk);
                coIter = coIter+1;
            elseif strcmp(D.Polarization{jj},'VV') || strcmp(D.Polarization{jj},'VH')
                sCx(:,cxIter) = D.Radar{jj,ii}(:,kk);
                cxIter = cxIter+1;
            end
            end
            D.Energy{jj,ii} = D.Radar{jj,ii}(:,kk).^2;
            iter = iter + 1;
        end
        
        for tt = floor(L/2):length(s)-ceil(L/2) %1:length(s)%
            c = s((tt-floor(L/2))+1:(tt+ceil(L/2)),:);
            % Compute Coherency
            c = sqrt(c.^2).*H;
            c1  = sum( (sum(c,2)).^2);  % Stacked Energy
            c2  = sum( sum(c.^2));      % Unstacked Energy 
            % Half Difference Between Out Energy and In Energy
            C(tt,kk) = 0.5.*(c1-c2);    % unnormalized cross-correlation sum
           if isCrossCoherence
            % Compute Cross-Coherency
            cCo = sCo((tt-floor(L/2))+1:(tt+ceil(L/2)),:);
            c = sqrt(cCo.^2).*H(:,[1,2]);
            c1  = sum( (sum(c,2)).^2);  % Stacked Energy
            c2  = sum( sum(c.^2));      % Unstacked Energy 
            % Half Difference Between Out Energy and In Energy
            Co(tt,kk) = 0.5.*(c1-c2);    % unnormalized cross-correlation sum
            
            % Compute Cross-Coherency
            cCx= sCx((tt-floor(L/2))+1:(tt+ceil(L/2)),:);

            c = sqrt(cCx.^2).*H(:,[1,2]);
            c1  = sum( (sum(c,2)).^2);  % Stacked Energy
            c2  = sum( sum(c.^2));      % Unstacked Energy 
            % Half Difference Between Out Energy and In Energy
            Cx(tt,kk) = 0.5.*(c1-c2);    % unnormalized cross-correlation sum
           end
            
        end
    end
    % Pad Image with Nearest Coherence
    padBot = ones(length(s)-tt,kk).*C(tt,:);
    padTop = ones(floor(L/2)-1,kk).*C(floor(L/2),:);
    C(1:floor(L/2)-1,:) = padTop;
    C(tt+1:end,:) = padBot;
    % Apply Smoothing and Normalize 
    CC = StakR(C,L);
    CC = CC./max(CC);
    D.Coherence{ii} = CC;
    
    if isCrossCoherence
    % Cross-Coherency HH-HV
    % Pad Image with Nearest Coherence
    padBot = ones(length(s)-tt,kk).*Co(tt,:);
    padTop = ones(floor(L/2)-1,kk).*Co(floor(L/2),:);
    Co(1:floor(L/2)-1,:) = padTop;
    Co(tt+1:end,:) = padBot;
    % Apply Smoothing and Normalize 
    CO = StakR(Co,L);
    CO = CO./max(CO);
    D.HHHVCoherence{ii} = CO;
    
    % Cross-Coherency VV-VH
    % Pad Image with Nearest Coherence
    padBot = ones(length(s)-tt,kk).*Cx(tt,:);
    padTop = ones(floor(L/2)-1,kk).*Cx(floor(L/2),:);
    Cx(1:floor(L/2)-1,:) = padTop;
    Cx(tt+1:end,:) = padBot;
    % Apply Smoothing and Normalize 
    CX = StakR(Cx,L);
    CX = CX./max(CX);
    D.VVVHCoherence{ii} = CX;
    end
    
    clear('s','sCo','sCx')
end
clear('L','H','C','CC','c','c1','c2','s','sCo','sCx','padTop','padBot',...
    'cCo','cCx','Cx','Cc','coIter','coPick','Co','cxIter')