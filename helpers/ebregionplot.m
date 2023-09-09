function h = ebregionplot(x, y, low, hi, colr, varargin)
p = inputParser;
addParameter(p,'ax', []);
parse(p,varargin{:});

%  h = ebregionplot(x,y,low,hi,plstr);
%
%  Plots gray error region on the current set of axes
%   (positioned underneath the last thing plotted)
%
%  Error region defined by [y-low, y+hi]
%
% Inputs: x, y  - ordinate and abscissa values
%   low = lower bound for error trace
%   hi = upper bound for error trace
%   colr = color for error region
%
% Outputs: h = handle to region;

if nargin < 5
    colr = .6*[1 1 1]; % default color
end

if isempty(p.Results.ax)
    ax = gca();
else
    ax = p.Results.ax;
end

% Reshape inputs into column vectors
x = x(:); y = y(:);  % make into column vector
low = low(:);
hi = hi(:);

% Query hold state
holdstate = ishold(ax);
hold on;

xx = [x;flipud(x)];
yy = [y-low; flipud(y+hi)];
fill(ax, xx,yy,colr,'edgecolor', colr,'facealpha',0.5,'edgealpha',0.);

% Place error surface under previous plots:
chldrn = get(ax, 'children');
if ~isempty(chldrn)
    set(ax, 'children', chldrn([2,1,3:end]));
end


if ~holdstate
    hold off;
end
