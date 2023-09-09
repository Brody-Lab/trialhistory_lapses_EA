function[] = plot_mean_RTs(data, varargin)


p = inputParser;

addParameter(p,'relative_to_mean', true);
addParameter(p, 'ploterrors', true);
addParameter(p, 'plotabs', true);

addParameter(p,'axHandle', []);
addParameter(p,'grid', 'off');
addParameter(p,'dataMarkerSize', 18);
addParameter(p,'dataLineWidth', 2);
addParameter(p, 'axisFontSize', 14);
addParameter(p, 'marker', 'none');
addParameter(p, 'dataLineStyle', '-');
addParameter(p, 'datatrans', false);

addParameter(p,'xlab', []);
addParameter(p,'ylab', 'Mean Reaction Time [s]');

parse(p,varargin{:});

if ~isempty(p.Results.axHandle)
    ax = p.Results.axHandle;
else
    ax = gca();
end

greens =  get_linear_cmap([46,139,87]./255, 10); % for correct trials
reds = get_linear_cmap([178,34,34]./255, 10);  % for error trials
if p.Results.datatrans
    col_hits = greens(8,:);
    col_err = reds(8,:);
else
    col_hits = greens(1,:);
    col_err = reds(1,:);
    
end

if p.Results.plotabs
    g = abs(data.gamma);
    xlim(ax, [0,5]);
    xticks = [0.5:1.5:4.5];
    xlab = '[Stimulus strength|';
else
    g = data.gamma;
    xlim(ax, [-5,5]);
    xticks = [-4.5:1.5:4.5];
    xlab = 'Stimulus strength';
end

if ~isempty(p.Results.xlab)
    xlab = p.Results.xlab;
end

xreg = unique(g);


for i = 1:length(xreg)
    idx_hits = g == xreg(i) & data.hits == 1;
    meanRT_hits(i) = mean(data.T(idx_hits));
    semRT_hits(i) = std(data.T(idx_hits))./sqrt(length(data.T(idx_hits)));
    
    idx_err = g == xreg(i) & data.hits == 0;
    meanRT_err(i) = mean(data.T(idx_err));
    semRT_err(i) = std(data.T(idx_err))./sqrt(length(data.T(idx_err)));
end

if p.Results.relative_to_mean
    meanRT_hits = meanRT_hits - 100;
    meanRT_err = meanRT_err - 100;
end

xreg = round(xreg,1);
errorbar(ax, xreg, meanRT_hits, semRT_hits, ...
    'color', col_hits,...
    'linewidth',p.Results.dataLineWidth,...
    'linestyle', p.Results.dataLineStyle, ...
    'marker', p.Results.marker,...
    'markersize',p.Results.dataMarkerSize,'capsize',0);

if p.Results.ploterrors
    errorbar(ax, xreg, meanRT_err, semRT_err, ...
        'color', col_err,...
        'linewidth', p.Results.dataLineWidth,...
        'linestyle', p.Results.dataLineStyle,...
        'marker', p.Results.marker,...
        'markersize', p.Results.dataMarkerSize,'capsize',0);
end


set(ax, 'XTick', xticks);
xtickangle(60)
grid(p.Results.grid)
xlabel(ax,xlab); 
ylabel(ax,p.Results.ylab);
ylim(ax, 'auto');

yticks(ax, 'auto');
box off;
set(ax, 'linewidth', 2)
set(ax,'TickDir','out')
set(ax, 'FontSize', p.Results.axisFontSize);
end