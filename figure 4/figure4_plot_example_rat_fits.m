master_startup;

example_rat = "X018";
saving = true;

datapath = fullfile(path.data_save,'pbupsRT/');
filenames =  get_reaction_time_filenames(datapath);
path_fitdata = fullfile(path.mat_save, 'saved_RTfits/');

%%

% load data and fit for example rat
rat_id = find(contains(filenames, example_rat));
display(filenames{rat_id})
load(datapath + filenames{rat_id});
avgdata.T = avgdata.T * 1000;  % convert to ms


fitdata = load(strcat(path_fitdata, filenames{rat_id}(1:end-4), "_goright.mat"));
fitdata.gamma = avgdata.gamma;


sim_idx = 1;
data.pokedR = fitdata.pred_choices(:,sim_idx)';
data.T = 1000*fitdata.pred_RTs(:,sim_idx)';
data.T(data.T>5000) = 5000;
data.gamma = avgdata.gamma;
data.hits = data.pokedR == sign(max(0, data.gamma));



%% Plot mean psychometric
figure('Name', example_rat);

ax_psy = subplot(1,4,1); hold on;
pbups_psych_gamma(avgdata, 'xreg', 'gamma',...
    'plotfit', false, ...
    'compute_fit', false,...
    'fitLineWidth', 1.5,...
    'axisFontSize', 14);

% fprintf("bootstrapping psychometric...\n")
% x = -4.8:0.1:4.8;
% [mean_y, ci_low, ci_up] = bootstrap_psych(bootstrap_psych_fit(data, 'nboots', 100), x);
% ebregionplot(x, mean_y, ci_low, ci_up, .4*[1 1 1], 'ax', ax_psy)

pbups_psych_gamma(data, 'xreg', 'gamma',...
    'plotfit', true, ...
    'compute_fit', true,...
    'plotdata', false,...
    'ploterrorbar', false,...
    'fitLineWidth', 1.5,...
    'axisFontSize', 14);
xticks([-4.5:3:4.5]);
xtickangle(0);

%% Plot conditioned psychometrics

ax_postc = subplot(1,4,2); hold on;
c_postright = get_linear_cmap([0 164 204]./255, 10);
c_postleft = get_linear_cmap([233 115 141]./255, 10);
x = -4.8:0.1:4.8;

plot_conditioned_psychometrics(avgdata, ax_postc, [],...
    'plot_fit', false,...
    'compute_fit', false,...
    'xreg', 'gamma', ...
    'axisFontSize', 14);


plot_conditioned_psychometrics(data, ax_postc, [],...
    'plot_fit', true,...
    'plot_data', false,...
    'compute_fit', true,...
    'xreg', 'gamma', ...
    'fitLineWidth', 1.5,...
    'axisFontSize', 14);

xticks([-4.5:3:4.5]);
xtickangle(0);

% load RT_bootstrap_cond_pych_params.mat;
% find the id for example rat
% rat_id = find(contains({rats(:).name}, example_rat));
% [mean_y, ci_low, ci_up] = bootstrap_psych(rats(rat_id).prm_fit.rc, x);
% ebregionplot(x, mean_y, ci_low, ci_up, c_postright(4,:), 'ax', ax_postc)
% [mean_y, ci_low, ci_up] = bootstrap_psych(rats(rat_id).prm_fit.lc, x);
% ebregionplot(x, mean_y, ci_low, ci_up, c_postleft(4,:), 'ax', ax_postc)

%% Plot correct and error RTs

ax_rt = subplot(1,4,3); cla; hold on;
plot_abs = false;

if plot_abs
    x = unique(abs(avgdata.gamma));
else
    x = unique(avgdata.gamma);
end

plot_mean_RTs(avgdata,...
    'dataLineStyle', '-',...
    'marker', 's',...
    'dataMarkerSize', 2.5,...
    'dataLineWidth', 1,...
    'plotabs', plot_abs, ...
    'relative_to_mean', false, ...
    'axHandle', ax_rt);

% okay now bootstrap RTs
% hits
greens =  get_linear_cmap([46,139,87]./255, 10); % for correct trials
[mean_y, ci] = bootstrap_mean_rts(fitdata, 1, 'plotabs', plot_abs);
ebregionplot(x, mean_y, ci(1,:), ci(2,:), greens(4,:), 'ax', ax_rt)
% errors
reds = get_linear_cmap([178,34,34]./255, 10);
[mean_y, ci] = bootstrap_mean_rts(fitdata, 0, 'plotabs', plot_abs);
ebregionplot(x, mean_y, ci(1,:), ci(2,:), reds(4,:), 'ax', ax_rt)

% hack to put data points back on top
plot_mean_RTs(avgdata,...
    'dataLineStyle', 'none',...
    'marker', 's',...
    'dataMarkerSize', 2.5,...
    'dataLineWidth', 1,...
    'plotabs', plot_abs, ...
    'relative_to_mean', false, ...
    'axHandle', ax_rt);
ylabel(ax_rt, 'Mean Reaction Time [ms]')
xlabel(ax_rt, 'Stimulus strength')
xticks([-4.5:3:4.5]);
xtickangle(0);
axis(ax_rt, 'square');




%% Plot post correct RTs

ax_postrt = subplot(1,4,4); cla; hold on;
x = unique(avgdata.gamma);

plot_bias_conditioned_RTs(avgdata, ...
    'axHandle', ax_postrt,...
    'dataLineStyle', '-',...
    'marker', 's', ...
    'dataLineWidth', 1,...
    'dataMarkerSize', 2.5);

% now bootstrap predictions
[mean_y, ci] = bootstrap_cond_rts(fitdata,...
    'do_hits', 1,...
    'do_right', 1);
ebregionplot(x, mean_y, ci(1,:), ci(2,:), c_postright(4,:), 'ax', ax_postrt)


[mean_y, ci] = bootstrap_cond_rts(fitdata,...
    'do_hits', 1,...
    'do_right', 0);
ebregionplot(x, mean_y, ci(1,:), ci(2,:), c_postleft(4,:), 'ax', ax_postrt)
xticks([-4.5:3:4.5]);
xtickangle(0);
axis(ax_postrt, 'square');


%% 
if saving
    set(gcf(), 'Units', 'inches', 'Position', [3,2,19,4]);
    savethisfig(gcf(), path.fig_save + "fig4_rtfits_initpt_" + example_rat)
end


%% bootstrap psychometric fits

function[mean_y, ci_low, ci_up] = bootstrap_psych(boot_psych, x)

nboots = length(boot_psych);
% okay make predictions for each gamma
pred_y = nan(nboots, length(x));
for i = 1:length(x)
    sigmoid = @(g0, g1, s, b, x) g0 + g1./(1+exp(-s.*(x - b)));
    pred_y(:,i) = sigmoid(boot_psych(:,1), boot_psych(:,2), boot_psych(:,3), boot_psych(:,4), x(i));
end
pred_y = sort(pred_y);

% 95 CI:
ci_idx = [floor(nboots*0.025), ceil(nboots*0.975)];
mean_y = mean(pred_y);
ci_low = mean_y - pred_y(ci_idx(1), :);
ci_up = pred_y(ci_idx(2), :) - mean_y;

end



%% bootstrap correct and error RT

% function[mean_y, ci] = bootstrap_mean_rts(data, do_hits, varargin)
%
% p = inputParser;
% addParameter(p, 'nboots', 1000);
% addParameter(p, 'plotabs', true);
% parse(p,varargin{:});
%
% nboots = p.Results.nboots;
% if p.Results.plotabs
%     g = abs(data.gamma);
% else
%     g = data.gamma;
% end
% xreg = unique(g);
%
% fprintf("bootstrapping RTs...\n")
% meanRT = nan(nboots, length(xreg));
% for i = 1:length(xreg)
%     idx = find(g == xreg(i) & data.hits == do_hits);
%     for n = 1:nboots
%         idx_boot = datasample(idx, length(idx));
%         meanRT(n,i) = mean(data.T(idx_boot));
%     end
% end
%
% % 95 CI:
% meanRT = sort(meanRT);
% ci_idx = [floor(nboots*0.025), ceil(nboots*0.975)];
% mean_y = mean(meanRT);
% ci = abs(meanRT(ci_idx, :) - mean_y);
%
% end


function[mean_y, ci] = bootstrap_mean_rts(data, do_hits, varargin)

p = inputParser;
addParameter(p, 'plotabs', true);
parse(p,varargin{:});

nboots = size(data.pred_choices,2);

if p.Results.plotabs
    g = abs(data.gamma);
else
    g = data.gamma;
end
xreg = unique(g);

fprintf("bootstrapping RTs...\n")
meanRT = nan(nboots, length(xreg));
for i = 1:length(xreg)
    for n = 1:nboots
        hits = data.pred_choices(:,n) == sign(max(0, data.gamma))';
        idx = g == xreg(i) & hits' == do_hits;
        meanRT(n,i) = mean(data.pred_RTs(idx,n)*1000);  % convert to ms
    end
end

% 95 CI:
meanRT = sort(meanRT);
ci_idx = [max(1,floor(nboots*0.025)), ceil(nboots*0.975)];
mean_y = mean(meanRT);
ci = abs(meanRT(ci_idx, :) - mean_y);

end


%% bootstrap post correct/error RTs
%
% function[mean_y, ci] = bootstrap_cond_rts(data, varargin)
%
% p = inputParser;
% addParameter(p, 'nboots', 1000);
% addParameter(p, 'do_hits', true);
% addParameter(p, 'do_right', true);
% parse(p,varargin{:});
%
% do_hits = p.Results.do_hits;
% do_right = p.Results.do_right;
% nboots = p.Results.nboots;
%
% % get post right choice trials
% pR = [0, (data.pokedR(1:end-1) == 1)];
% pH = [0, (data.hits(1:end-1) == 1)];
% g = data.gamma;
% xreg = unique(g);
%
% fprintf("bootstrapping RTs...\n")
% meanRT = nan(nboots, length(xreg));
% for i = 1:length(xreg)
%     idx = find(g == xreg(i) & pR == do_right & pH == do_hits & data.hits == 1);
%     for n = 1:nboots
%         idx_boot = datasample(idx, length(idx));
%         meanRT(n,i) = mean(data.T(idx_boot));
%     end
% end
%
% meanRT = meanRT - mean(mean(meanRT));
%
% % 95 CI:
% meanRT = sort(meanRT);
% ci_idx = [floor(nboots*0.025), ceil(nboots*0.975)];
% mean_y = mean(meanRT);
% ci = abs(meanRT(ci_idx, :) - mean_y);
%
% end

function[mean_y, ci] = bootstrap_cond_rts(data, varargin)

p = inputParser;
addParameter(p, 'do_hits', true);
addParameter(p, 'do_right', true);
parse(p,varargin{:});

do_hits = p.Results.do_hits;
do_right = p.Results.do_right;

nboots = size(data.pred_choices,2);

% get post right choice trials

g = data.gamma;
xreg = unique(g);

fprintf("bootstrapping RTs...\n")
meanRT = nan(nboots, length(xreg));
for i = 1:length(xreg)
    for n = 1:nboots
        pR = [0, (data.pred_choices(1:end-1,n) == 1)'];
        hits = data.pred_choices(:,n)' == sign(max(0, data.gamma));
        pH = [0, (hits(1:end-1) == 1)];
        idx = find(g == xreg(i) & pR == do_right & pH == do_hits & hits == 1);
        meanRT(n,i) = mean(data.pred_RTs(idx, n)*1000);
    end
end

meanRT = meanRT - mean(mean(meanRT));

% 95 CI:
meanRT = sort(meanRT);
ci_idx = [max(1,floor(nboots*0.025)), ceil(nboots*0.975)];
mean_y = mean(meanRT);
ci = abs(meanRT(ci_idx, :) - mean_y);

end
