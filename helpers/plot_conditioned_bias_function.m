function[] = plot_conditioned_bias_function(data, varargin)


p = inputParser;

addParameter(p,'nquantiles', 6);
addParameter(p,'axHandle', []);
addParameter(p,'grid', 'off');
addParameter(p,'dataMarkerSize', 18);
addParameter(p,'dataLineWidth', 2);
addParameter(p, 'axisFontSize', 14);
addParameter(p, 'datatrans', false);
addParameter(p,'xlab', 'Reaction time [s]');
addParameter(p,'ylab', 'Choice repetition bias');

parse(p,varargin{:});

if ~isempty(p.Results.axHandle)
    ax = p.Results.axHandle; hold on;
else
    ax = gca(); hold on;
end

if p.Results.datatrans
    col = [0.6, 0.6, 0.6];
else
    col = 'k';
end


% post correct trials
idx = [false, data.hits(1:end-1) == 1];
idx = idx & data.hits == 1;

% was choice a repetition
rep = [0, data.pokedR(1:end-1) == data.pokedR(2:end)];

% divide RTs into quantiles
rtquantiles = [min(data.RT(idx))-eps(), quantile(data.RT(idx), p.Results.nquantiles), max(data.RT(idx)+eps())];
rtcenters = (rtquantiles(1:end-1) + rtquantiles(2:end))/2;

% compute mean repetation bias for that RT quantile
bins = discretize(data.RT(idx), rtquantiles);
vals = grpstats(rep(idx), bins);

% plot
plot(ax, rtcenters, vals,...
    'color', col,...
    'linewidth', p.Results.dataLineWidth);
scatter(ax, rtcenters, vals, p.Results.dataMarkerSize, 'filled',...
    'MarkerFacecolor', col)
hline(0.5, 'k:')



% format plot
grid(p.Results.grid)
xlabel(ax,p.Results.xlab); ylabel(ax,p.Results.ylab);
ylim(ax, 'auto');
% xlim(ax, [0, 1.0]);
yticks(ax, 'auto');
box off;
set(ax, 'linewidth', 2)
set(ax,'TickDir','out')
set(ax, 'FontSize', p.Results.axisFontSize);
end
