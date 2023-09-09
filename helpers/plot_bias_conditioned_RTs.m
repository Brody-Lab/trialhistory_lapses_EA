function[] = plot_bias_conditioned_RTs(data, varargin)


p = inputParser;

addParameter(p,'nquantiles', 6);
addParameter(p,'axHandle', []);
addParameter(p,'grid', 'off');
addParameter(p,'dataMarkerSize', 18);
addParameter(p,'dataLineWidth', 2);
addParameter(p, 'marker', 'none');
addParameter(p, 'dataLineStyle', '-');
addParameter(p, 'axisFontSize', 14);
addParameter(p, 'datatrans', false);
addParameter(p, 'plot_uncond', false);
addParameter(p,'xlab', 'Stimulus strength');
addParameter(p,'ylab', '% change from mean RT');

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

c_postright = get_linear_cmap([0 164 204]./255, 10);
c_postleft = get_linear_cmap([233 115 141]./255, 10);

if p.Results.datatrans
    col_R = c_postright(7,:);
    col_L = c_postleft(7,:);
else
    col_R = c_postright(1,:);
    col_L = c_postleft(1,:);
end

% get post right choice trials
pR = [0, (data.pokedR(1:end-1) == 1)];
pH = [0, (data.hits(1:end-1) == 1)];
g = data.gamma;
xreg = unique(g);

for i = 1:length(xreg)
    idx = g == xreg(i) & pR == 1 & pH == 1 & data.hits == 1;
    meanRT_post_right(i) = mean(data.T(idx));
    semRT_post_right(i) = std(data.T(idx))./sqrt(length(data.T(idx)));
    
    idx = g == xreg(i) & pR == 0 & pH == 1 & data.hits == 1;
    meanRT_post_left(i) = mean(data.T(idx));
    semRT_post_left(i) = std(data.T(idx))./sqrt(length(data.T(idx)));
end

mean_resp_R = meanRT_post_right - mean(meanRT_post_right);
mean_resp_L = meanRT_post_left - mean(meanRT_post_left);

errorbar(ax, xreg(xreg<0), mean_resp_R(xreg<0), semRT_post_right(xreg<0), ...
    'color', col_R,...
    'linewidth',p.Results.dataLineWidth,...
    'linestyle', p.Results.dataLineStyle,...
    'marker', p.Results.marker,...
    'markersize',p.Results.dataMarkerSize,'capsize',0);

errorbar(ax, xreg(xreg>0), mean_resp_R(xreg>0), semRT_post_right(xreg>0), ...
    'color', col_R,...
    'linewidth',p.Results.dataLineWidth,...
    'linestyle', p.Results.dataLineStyle,...
    'marker', p.Results.marker,...
    'markersize',p.Results.dataMarkerSize,'capsize',0);


errorbar(ax, xreg(xreg<0), mean_resp_L(xreg<0), semRT_post_right(xreg<0), ...
    'color', col_L,...
    'linewidth',p.Results.dataLineWidth,...
    'linestyle', p.Results.dataLineStyle,...
    'marker', p.Results.marker,...
    'markersize',p.Results.dataMarkerSize,'capsize',0);

errorbar(ax, xreg(xreg>0), mean_resp_L(xreg>0), semRT_post_right(xreg>0), ...
    'color', col_L,...
    'linewidth',p.Results.dataLineWidth,...
    'linestyle', p.Results.dataLineStyle,...
    'marker', p.Results.marker,...
    'markersize',p.Results.dataMarkerSize,'capsize',0);

if p.Results.plot_uncond
    
    for i = 1:length(xreg)
        idx = g == xreg(i) & pH == 1 & data.hits == 1;
        meanRT(i) = mean(data.T(idx));
        semRT(i) = std(data.T(idx))./sqrt(length(data.T(idx)));
    end
    meanRT = meanRT - mean(meanRT);
    
    errorbar(ax, xreg, meanRT, semRT, ...
        'color', col, 'linestyle', ':', ...
        'linewidth',p.Results.dataLineWidth,...
        'linestyle', p.Results.dataLineStyle,...
        'marker', p.Results.marker,...
        'markersize',p.Results.dataMarkerSize,'capsize',0);
end

grid(p.Results.grid)
xlabel(ax,p.Results.xlab); ylabel(ax,p.Results.ylab);
ylim(ax, 'auto');
xlim(ax, [-5, 5]);
yticks(ax, 'auto');
box off;
set(ax, 'linewidth', 2)
set(ax,'TickDir','out')
set(ax, 'FontSize', p.Results.axisFontSize);

end
