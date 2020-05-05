%% CMP Velocity Analysis
% Choose Polarization
polAns = questdlg('Choose the CMP Polarization.','CMP Polarization','HH','HV','HH');
MD.polIx = find(strcmp({polAns},MD.Polarization));
% Choose Number of Horzons to Pick
prompt = {'Number of Direct Wave Horizons','Numer of Reflection Horizons','Direct Wave Norm', 'Reflection Norm' };
dlgtitle = 'Input Horizons';
dfaultinput = {'1','1','L2','L2'};
dim = 1;
horizons = inputdlg(prompt,dlgtitle,dim,dfaultinput);
MD.ndir = str2num(horizons{1}); % Number of Direct Horizons
MD.nref = str2num(horizons{2}); % Number of Reflection Horizons
% Choose Inversion norm
MD.dirNorm = horizons{3};
MD.refNorm = horizons{4};

try close(998)
    close(998)
catch
end
% Pick Direct Wave
for kk = 1:MD.ndir
    picks = ManualPicker(D.CMP{MD.polIx},'zero',3,MD.dt,MD.offsetArray,MD.TimeAxis);
    D.dirTpick{kk} = picks;
end

% Pick Ground Reflection
for kk = 1:MD.nref
    picks = ManualPicker(D.CMP{MD.polIx},'zero',3,MD.dt,MD.offsetArray,MD.TimeAxis);
    D.refTpick{kk} = picks;
end

clear('dim','dfaultinput','prompt','polAns','horizons','picks')
pause(5)
close(998)
%% BootStrapping Velocity Analysis
MCix = 1:length(MD.offsetArray);
nMC = 250; % Number of MC Simulations
for kk = 1:nMC
    binIx = datasample(MCix,round(0.5.*length(MD.offsetArray)));
    % Least Squares Velocity Analysis
    % Direct Wave Analysis
    for hh = 1:MD.ndir
        G = [ones(length(MD.offsetArray(binIx)),1),MD.offsetArray(binIx)'];
        d = D.dirTpick{hh}(binIx);
        if strcmp(MD.dirNorm,'L2')
            % L2 norm
            m = G\d;
        else
            %L1 norm
            p = 1; % p-Norm
            tolr = 1e-2; % Residual Tolerance
            tolx = 5e-4; % Convergence Tolerance
            maxiter = 100; % Maximm Iterations
            m = irls(G,d,tolr,tolx,p,maxiter);
        end
        % Save Intercept Time, Velocity, Wavelength, and Density
        D.dirT0{hh}(kk) = m(1);
        D.dirV{hh}(kk) = 1/m(2);
        D.dirZ{hh}(kk) = D.dirV{hh}(kk)./MD.f0GHz;
        D.dirRho{hh}(kk) = DryCrim(D.dirV{hh}(kk));
    end
    % Reflecttion Velocity Analysis
    for hh = 1:MD.nref
        G = [ones(length(binIx),1),(MD.offsetArray(binIx))'.^2];
        d = [(D.refTpick{hh}(binIx) - D.dirT0{1}(kk)).^2];
        if strcmp(MD.refNorm,'L2')
            % L2 norm
            m = G\d;
        else
            % L1 norm
            p = 1; % p-Norm
            tolr = 1e-2; % Residual Tolerance
            tolx = 5e-4; % Convergence Tolerance
            maxiter = 100; % Maximm Iterations
            m = irls(G,d,tolr,tolx,p,maxiter);
        end
        % Save Intercept Time, Velocity, Depth, and Density
        D.refT0{hh}(kk) = sqrt(m(1));
        D.refV{hh}(kk) = 1./sqrt(m(2));
        D.refZ{hh}(kk) = D.refT0{hh}(kk).*D.refV{hh}(kk)./2;
        D.refRho{hh}(kk) = DryCrim(D.refV{hh}(kk));
    end
end
clear('G','d','m','ix','picks','binIx')