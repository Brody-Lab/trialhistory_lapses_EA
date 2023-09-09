
master_startup;
saving = true;


datapath = fullfile(path.data_save,'pbupsRT/');
filenames =  get_reaction_time_filenames(datapath);
[gamma, pokedR, hits] = deal([]);


%% Initialize choice trial-history figure

f = figure(); hold on;
ax_psych = subplot(1,8,[1,2]); hold on;
ax_postcorr = subplot(1,8,[3,4]); hold on;
for ax = 1:4
    ax_prm(ax) = subplot(1,8, 4+ ax); hold on;
    ax_prm(ax).PlotBoxAspectRatio =  [0.275 0.75 0.75];
end

%% loop thorugh rats and plot psychometrics (history conditioned and not)

for i = 1:length(filenames)

    load(datapath +  filenames{i});
    
    % Psychometrics
    pbups_psych_gamma(avgdata, 'axHandle', ax_psych, ...
        'fitLineWidth', 1.5, ...
        'fitLineColor', [0.8 0.8 0.8], ...
        'dataLineColor', [0.8 0.8 0.8], ...
        'errorbarColor',[0.8 0.8 0.8], ...
        'dataFaceColor', [0.8, 0.8, 0.8]);
    
    % conditioned psychometrics
    prm{i} = plot_conditioned_psychometrics(avgdata, ax_postcorr, [], 'trans',true);

    % reformat psychometric parameters for plotting param scatters
    rats(i).name  = avgdata.ratname{1};
    rats(i).prm.rc = unpack_psychparams(prm{i}, 'beta_rc');
    rats(i).prm.lc = unpack_psychparams(prm{i}, 'beta_lc');
    rats(i).ntrials = length(avgdata.hits);
    rats(i).accuracy = sum(avgdata.hits)/rats(i).ntrials;
    
    
    % concatenate for meta rat plots
    gamma = [gamma avgdata.gamma];
    pokedR = [pokedR avgdata.pokedR];
    hits = [hits avgdata.hits];
   
end

%% plot meta rat data

% settle gammas across rats
targets = -4.5:4.5;
targets = setdiff(setdiff(targets, 2.5), -2.5);
gamma = interp1(targets, targets, gamma, 'nearest');

data = [];
data.pokedR = pokedR;
data.gamma = gamma;
data.hits = hits;

% mean psychometric
pbups_psych_gamma(data, ...
    'fitLineWidth', 2.5,...
    'fitLineColor', [0 0 0],...
    'plotdata', false, 'ploterrorbar', false,...
    'axHandle', ax_psych);

% conditioned psychometrics
plot_conditioned_psychometrics(data, ax_postcorr, [], 'plot_data', false);
ylabel(ax_postcorr, [])
yticks(ax_postcorr, [])



%% plot parameter scatters

c_postright = get_linear_cmap([0 164 204]./255, 10);
c_postleft = get_linear_cmap([233 115 141]./255, 10);
c_pR = c_postright(1,:);
c_pL = c_postleft(1,:);

d = reorganize_prm(rats, 'prm');
params = {'sens', 'bias', 'left_lapse', 'right_lapse'};
titles = {'SENSITIVITY', 'BIAS', 'LEFT LAPSE', 'RIGHT LAPSE'};

fprintf('\n ==== \n')
for i = 1:length(params)
    
    plotdata = [d.rc.(params{i})(1:end)' d.lc.(params{i})(1:end)'];
    
    % center the bias parameter
    if i == 2
        plotdata = plotdata - mean(plotdata,2);
        plotdata
    end
    
    plot_param_scatters(plotdata, ax_prm(i),  c_pR, c_pL);
    title(ax_prm(i), titles{i});
    
    % plotting mean again, so it is on top
    plot(ax_prm(i), [1 1.5], median(plotdata), '-o',...
        'markerfacecolor', 'k',...
        'markersize',10,...
        'LineWidth', 2.5, ...
        'color',  'k');
    scatter(ax_prm(i), 1, median(plotdata(:,1)), 50, c_pR, 'filled');
    scatter(ax_prm(i), 1.5, median(plotdata(:,2)), 50, c_pL, 'filled');
    
    [p,h] = ranksum(plotdata(:,1), plotdata(:,2));
    if h == 0
        fprintf('\n %s is NOT significantly different with p-value %d\n', titles{i}, p)
    else
        fprintf('\n %s IS significantly different with p-value %d\n', titles{i}, p)
    end
    
end

ylim(ax_prm(1), [0.9, 1.5])
ylim(ax_prm(2), [-0.7, 0.7])
ylim(ax_prm(3), [0, 0.35])
ylim(ax_prm(4), [0, 0.35])
drawnow;


%%
if saving
    set(0,'DefaultFigureWindowStyle', 'normal')
    set(f, 'Units', 'inches', 'Position', [3,2,19,8])
    savethisfig(f, 'fig4_rtchoicehistoryeffects')
end



