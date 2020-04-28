function [pout] = distributeDensity(E,N,rho,X,Y,r,p)
%distributeDensity used distance weighting to calculate the average density
%along the GPR transect from the nearest pit locations (X,Y)
if nargin<7
    p = 2;
end
if nargin < 6
    r = 1000;
end
dist = zeros(length(E),length(X));
for kk = 1:length(E)
dist(kk,:) = sqrt((E(kk)-X).^2+(N(kk)-Y).^2);
end
isDist = dist<r;
% Inverse Distance Weighting
pout = zeros(length(dist),1);
for kk = 1:length(dist)
    pitIx = isDist(:,kk);
    w = 1./dist(pitIx,kk).^p;
    pout(kk) = sum(w.*rho(pitIx))./sum(w);
end
end
