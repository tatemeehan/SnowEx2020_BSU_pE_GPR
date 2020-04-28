function picks = ManualPicker(varargin)
% 
% INPUT:     varargin{1} signal = matrix of input data
%            varargin{2} method = 'min', 'max', or 'zero' for peak snap
%            varargin{3} window = [samples] number of samples to look 
%                                   +- the raw pick for a peak to snap to
%            varargin{4} dt     = [ns] sample rate
%            varargin{5} x      = [m] x axis (Distance or Offset)
%            varargin{6} y      = [ns] y axis (Travel-time)
% 
% OUTPUT:     picks = [index] interpolated, rounded pick values

signal      = varargin{1};
ntr         = length(signal(1,:)); % number of traces

if length(varargin) == 3
    window  = varargin{3};
    method  = varargin{2};
    dt      = .1;
    x       = 1:ntr;
    y       = 1:size(signal,1);
    
elseif length(varargin) == 2
    method  = varargin{2};
    window  = 3;
    dt      = .1;
    x       = 1:ntr;
    y       = 1:size(signal,1);
elseif length(varargin) == 1
    window  = 3;
    method  = 'zero';
    dt      = .1;
    x       = 1:ntr;
    y       = 1:size(signal,1);
elseif length(varargin) == 4
    window  = varargin{3};
    method  = varargin{2};
    dt      = varargin{4};
    x       = 1:ntr;
    y       = 1:size(signal,1);
elseif length(varargin) == 5
    window  = varargin{3};
    method  = varargin{2};
    dt      = varargin{4};
    x       = varargin{5};
    y       = 1:size(signal,1);
elseif length(varargin) == 6
    window  = varargin{3};
    method  = varargin{2};
    dt      = varargin{4};
    x       = varargin{5};
    y       = varargin{6};
end


% Plot the data
% if any(findobj(allchild(0), 'flat', 'type', 'figure') == h)
figure(998); 
c = get(gca,'Children');%clf
try c(length(c)).CData;
catch
% imagesc(normalize(signal)); hold on
imagesc(x,y,signal);colormap(bone); hold on
end
% Pick the surface
% ginput rules: left mouse click = pick, any key = quit
button = 1;
kk = 1;
while button == 1
[tmpx, tmpy, button] = ginput(1);
if ~isempty(tmpx)
rpx(kk) = tmpx; rpy(kk) = tmpy;
plot(rpx(kk),rpy(kk),'rx')
kk = kk+1;
end
end
rpx=rpx(:);rpy = rpy(:);
nrp = length(rpx);
rp = [rpx,rpy];
% Sort as ascending x values
[sorted,IxS]    = sort(rp(:,1),'ascend');
RawPicks(:,1)   = sorted;
RawPicks(:,2)   = rp(IxS,2);

% Now get only unique values 
[RawX,IxU]      = unique(RawPicks(:,1));
RawY            = RawPicks(IxU,2);

% Fill in
RawPicksAll     = interp1(RawX,RawY,x,'pchip');
IxN             = find(isnan(RawPicksAll));
IxT             = diff(IxN)-1;
mid             = find(IxT>0);
RawPicksAll(1:IxN(mid))     = RawPicksAll(mid+1);
RawPicksAll(IxN(mid+1):end) = RawPicksAll(IxN(mid+1)-1);

% Get the actual peaks within the window
picks = zeros(ntr,1);
for n = 1:ntr
    if nargin >= 5
        ixRange        = int8((round(RawPicksAll(n),1)-window.*dt)./dt):...
            int8((round(RawPicksAll(n),1)+window.*dt)./dt);
    else
        ixRange        = int8(round(RawPicksAll(n))-window):...
            int8(round(RawPicksAll(n))+window);
    end
    % Trave-Time Picks Snapped to Method
    [picks(n)]	= snapPicks(signal(:,n),RawPicksAll(n),ixRange,dt,window,method);
end
% Refresh Axis Children
c = get(gca,'Children');
% Remove Raw Picks from figure(998)
delete(c(1:nrp))
% Refresh Axis Children
c = get(gca,'Children');
% Plot Interpolated Picks
plot(x,picks,'rx')

end