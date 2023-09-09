
master_startup;
saving = true;

datapath = fullfile(path.data_save,'pbupsRT/');
filenames =  get_reaction_time_filenames(datapath);
[gamma, pokedR, T, RT, hits] = deal([]);


nquantiles = 6;


%% Initialize figure for plotting

fig = figure('Name', 'data'); hold on;
ax_rts = subplot(1,6,1); hold on;
ax_pc = subplot(1,6,[2,3]); hold on;
ax_crf = subplot(1,6,[4,5,6]); hold on;

%%
figure();
ax_rt_mean = gca(); hold on;
% ax_rt_mean.PlotBoxAspectRatio =  [0.275 0.7 0.75];

% loop through rats and plot reaction times with various conditioning
for k = 1:length(filenames)
        
    load(datapath + filenames{k});
    
    data = [];
    data.pokedR = avgdata.pokedR;
    data.gamma = avgdata.gamma;
    data.hits = avgdata.hits;
    data.RT = avgdata.T;
    
    % demean wrt side-specific non-decision times
    data.T(data.pokedR == 1) = avgdata.T(data.pokedR == 1) /mean(avgdata.T(data.pokedR == 1));
    data.T(data.pokedR == 0) = avgdata.T(data.pokedR == 0) /mean(avgdata.T(data.pokedR == 0));
    data.T = 100.*data.T;
    
    % concatenate for meta rat plots
    gamma = [gamma avgdata.gamma];
    pokedR = [pokedR avgdata.pokedR];
    hits = [hits avgdata.hits];
    T = [T data.T];
    RT = [RT data.RT];
    
    % Fraction change in RTs
    plot_mean_RTs(data, 'axHandle', ax_rts,...
        'plotabs', true,...
        'datatrans', true, ...
        'dataLineWidth', 1.5,...
        'ploterrors', true,...
        'ylab', '% change from mean RT');
    
    % MEAN RTs
    scatter(ax_rt_mean, 1, mean(avgdata.T), 50, 'filled','markerFaceColor', 'k', 'MarkerEdgeColor',[1 1 1], 'LineWidth',1.)
    errorbar(ax_rt_mean, 1, mean(avgdata.T), std(avgdata.T)/sqrt(length(avgdata.T)), ...
        'LineStyle', 'none',...
        'color', 'k',...
        'markerfacecolor','k',...
        'linewidth',1.5,...
        'markersize',80,'capsize',0);
    
    % (initial point/choice) conditioned RTs
    plot_bias_conditioned_RTs(data, 'axHandle', ax_pc, ...
        'datatrans', true, 'dataLineWidth', 1);
    
    
    % conditional bias functions
    plot_conditioned_bias_function(data, 'axHandle', ax_crf,...
        'datatrans', true, ...
        'dataLineWidth', 1, ...
        'nquantiles', nquantiles,...
        'dataMarkerSize', 12);
    
    
    
end

% format mean RT plot
ylabel(ax_rt_mean, 'Mean RT [s]');
ylim(ax_rt_mean, [0.125, 0.3]);
yticks(ax_rt_mean, 'auto');
xticks(ax_rt_mean, []);
set(ax_rt_mean, 'linewidth', 2.5)
set(ax_rt_mean,'TickDir','out')
set(ax_rt_mean, 'FontSize', 18);
box off;
if saving
    set(gcf(), 'Units', 'inches', 'Position', [3,2,5,4]);
    savethisfig(gcf(), path.fig_save + 'fig4_meanrt_data')
    close(gcf);
end

%% plot meta rat

% settle gammas across rats
targets = -4.5:4.5;
targets = setdiff(setdiff(targets, 2.5), -2.5);
gamma = interp1(targets, targets, gamma, 'nearest');

data.pokedR = pokedR;
data.gamma = gamma;
data.hits = hits;
data.T = T;
data.RT = RT;

% mean fractional chane in RTs across rats
plot_mean_RTs(data, 'axHandle', ax_rts,...
    'ploterrors', true,...
    'plotabs', true,...
    'datatrans', false, ...
    'dataLineWidth', 2.5,...
    'ylab', '% change from mean RT');
ylim(ax_rts, [-25, 12]);

plot_bias_conditioned_RTs(data, 'axHandle', ax_pc, ...
    'datatrans', false, 'plot_uncond', true);
ylim(ax_pc, [-25 25]);

plot_conditioned_bias_function(data, 'axHandle', ax_crf,...
    'datatrans', false, ...
    'dataLineWidth', 2, ...
    'nquantiles', nquantiles, ...
    'dataMarkerSize', 20)
ylim(ax_crf, [0.4, 0.8])
plot(ax_crf, xlim(ax_crf), [0.5, 0.5], 'k:');

if saving
    set(fig, 'Units', 'inches', 'Position', [3,2,19,4]);
    savethisfig(fig, path.fig_save + 'fig4_rtdata')
end


%% Set some common parameters for simulation plots

bound = 4.9;  % 4.9
sigma_sens = 2.5; % 2.5
ntrials = 120000;

%%

% Initial point
pc = 1.0;   % 1.0
tc = 0.7;   % 0.7
pe = 0.2;   % 0.2
te = 0.1;   % 0.1
plot_simulation_rt_bias('history_params', make_history_params(pc, pe, tc, te),...
    'stim_type', 'discrete',...
    'hist_effect', 'initial_pt',...
    'history_biased', true,...
    'ntrials', ntrials,...
    'bound', bound,...
    'seed', 1,...  %1
    'sigma_sens', sigma_sens,...
    'nquantiles', nquantiles);

if saving
    set(gcf(), 'Units', 'inches', 'Position', [3,2,19,4]);
    savethisfig(gcf(), path.fig_save + 'fig4_rt_initpt')
end

%%

% No history
[pc, pe, tc, te] = deal(0);
plot_simulation_rt_bias('history_params', make_history_params(pc, pe, tc, te),...
    'stim_type', 'discrete',...
    'hist_effect', 'initial_pt',...
    'history_biased', false,...
    'ntrials', ntrials,...
    'bound', bound,...
    'seed', 2,...
    'sigma_sens', sigma_sens,...
    'nquantiles', nquantiles);

if saving
    set(gcf(), 'Units', 'inches', 'Position', [3,2,19,4]);
    savethisfig(gcf(), path.fig_save + 'fig4_rt_nohistory')
end


%%


function[] = plot_simulation_rt_bias(varargin)


p = inputParser;
addParameter(p, "stim_type", "discrete", @(x) ismember(x, {'continuous', 'discrete'}));
addParameter(p, "hist_effect", "initial_pt", @(x) ismember(x, {'initial_pt', 'drift', 'both'}));
addParameter(p, "nquantiles", 6);
addParameter(p, "ntrials", 80000);
addParameter(p, "seed", 1);
addParameter(p, "bound", 4);
addParameter(p, "sigma_sens", 2.0);
addParameter(p, "mus", sort([-4.5:1:4.5]));
addParameter(p, "history_biased", true);
addParameter(p, "history_params", make_history_params(0.25, -0.25, 0.8, 0.8));
parse(p, varargin{:});

% make new parameters based on inputs
params = fields(p.Results);
for i = 1:length(params)
    eval([params{i} '=p.Results.' params{i} ';' ])
end

if ~history_biased 
    figure('Name', 'No history');
else
    figure('Name', hist_effect);
end

data = simDDM('stim_type', stim_type, ...
    'bound', bound,...
    'ntrials', ntrials,...
    'sigma_sens', sigma_sens,...
    'hist_effect', hist_effect,...
    'history_biased', true,...
    'history_params', history_params,...
    'mus', mus,...
    'seed', seed);

data.RT = data.T + random('InverseGaussian',0.01,2,ntrials,1)';
data.T =100.*data.RT/mean(data.RT);

ax1 = subplot(1,6,1); hold on;
plot_mean_RTs(data, 'axHandle', ax1,...
    'ploterrors', true,...
    'plotabs', true,...
    'datatrans', false, ...
    'dataLineWidth', 2.5,...
    'ylab', '% change from mean RT');
ylim([-20, 12]);

ax2 = subplot(1,6,[2,3]); hold on;
plot_bias_conditioned_RTs(data, 'axHandle', ax2, ...
    'datatrans', false, 'plot_uncond', true);
ylim([-25 25]);

ax3 = subplot(1,6,[4,5,6]); hold on;
plot_conditioned_bias_function(data, 'axHandle', ax3,...
    'datatrans', false, ...
    'dataLineWidth', 2, ...
    'nquantiles', nquantiles)
ylim([0.4, 0.8])

% ax4 = subplot(1,4,4); hold on;
% plot_conditioned_psychometrics(data, gca(), [], 'plot_data', false);
end