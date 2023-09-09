%%
master_startup;

datapath = fullfile(path.data_save,'pbupsRT/');
filenames =  get_reaction_time_filenames(datapath);

saving = 1;
example_rat = 'X018';

%%
load(path.mat_save + "RT_bootstrap_cond_pych_params.mat");

[ci_bias, ci_lapse, ci_fitbias, ci_fitlapse] = deal(nan(length(filenames),2));
[mean_bias, mean_lapse, mean_fitbias, mean_fitlapse] = deal(nan(length(filenames),1));


for i = 1:6
    
    display(filenames{i});
    
    % compute 95% CIs
    [ci_bias(i,:), ci_lapse(i,:), mean_bias(i), mean_lapse(i)] = compute_lapsebias_CI(rats(i).prm);
    [ci_fitbias(i,:), ci_fitlapse(i,:), mean_fitbias(i), mean_fitlapse(i)] = compute_lapsebias_CI(rats(i).prm_fit);
    
  
    
end

example_id = find(contains(filenames, example_rat));

%%
figure();

col = [43, 123, 186]./255;
ax1 = subplot(1,2,1); hold on;
errorbar(mean_bias, mean_fitbias,...
    ci_fitbias(:,1),...
    ci_fitbias(:,2),...
    ci_bias(:,1),...
    ci_bias(:,2),...
    'd',...
    'capsize',0,...
    'MarkerFaceColor', col,...
    'MarkerEdgeColor', col,...
    'Color', col,...
    'linewidth',1)
scatter(mean_bias(example_id), mean_fitbias(example_id), '*');
format_prm_scatter(ax1, -1.5, 0.5);

ax2 = subplot(1,2,2); hold on;
errorbar(mean_lapse, mean_fitlapse,...
    ci_fitlapse(:,1),...
    ci_fitlapse(:,2),...
    ci_lapse(:,1),...
    ci_lapse(:,2),...
    'd',...,
    'capsize',0,...
    'MarkerFaceColor', col,...
    'MarkerEdgeColor', col,...
    'Color', col,...
    'linewidth',1)
scatter(mean_lapse(example_id), mean_fitlapse(example_id), '*');
format_prm_scatter(ax2, -0.05, 0.6);



%%
load(path.mat_save + "RT_bootstrap_cond_pych_params_lapsemod.mat");

[ci_bias, ci_lapse, ci_fitbias, ci_fitlapse] = deal(nan(length(filenames),2));
[mean_bias, mean_lapse, mean_fitbias, mean_fitlapse] = deal(nan(length(filenames),1));


for i = 1:6
    
    display(filenames{i});
    
    % compute 95% CIs
    [ci_bias(i,:), ci_lapse(i,:), mean_bias(i), mean_lapse(i)] = compute_lapsebias_CI(rats(i).prm);
    [ci_fitbias(i,:), ci_fitlapse(i,:), mean_fitbias(i), mean_fitlapse(i)] = compute_lapsebias_CI(rats(i).prm_fit);
    
end

%%

% figure()
col = [0.6, 0.6, 0.6];
ax1 = subplot(1,2,1); hold on;
errorbar(mean_bias, mean_fitbias,...
    ci_fitbias(:,1),...
    ci_fitbias(:,2),...
    ci_bias(:,1),...
    ci_bias(:,2),...
    's',...
    'capsize',0,...
    'MarkerFaceColor', col,...
    'MarkerEdgeColor', col,...
    'Color', col,...
    'linewidth',1)
format_prm_scatter(ax1, -1.5, 0.2);

ax2 = subplot(1,2,2); hold on;
errorbar(mean_lapse, mean_fitlapse,...
    ci_fitlapse(:,1),...
    ci_fitlapse(:,2),...
    ci_lapse(:,1),...
    ci_lapse(:,2),...
    'sc',...
    'capsize',0,...
    'MarkerFaceColor', col,...
    'MarkerEdgeColor', col,...
    'Color', col,...
    'linewidth',1)
format_prm_scatter(ax2, -0.05, 0.5);

if saving
    set(gcf(), 'Units', 'inches', 'Position', [3,2,10,4]);
    savethisfig(gcf(), path.fig_save + "fig4_lapsebias")
end

%% MODEL COMPARISON

filenames =  get_reaction_time_filenames(fullfile(path.data_save,'pbupsRT/'));
path_fitdata = fullfile(path.mat_save, 'saved_RTfits/');


% bic = label_nparams*log(label_ntrials) - 2*label_ll;
for i = 1:length(filenames)
    
    display(filenames{i})
%     data = load([path_pkgdata, filenames{i}]);
    fitdata = load(strcat(path_fitdata, filenames{i}(1:end-4), "_lapsemod_goright.mat"));
    
    
    lapsemod.loglike(i) = fitdata.loglik;
    lapsemod.fit_params(i) = sum(fitdata.fit);
    lapsemod.ntrials(i) = length(fitdata.initial_pt);
    lapsemod.bic(i) = sum(fitdata.fit)*log(length(fitdata.initial_pt)) - 2*fitdata.loglik;

end


for i = 1:length(filenames)
    
    display(filenames{i})
%     data = load([path_pkgdata, filenames{i}]);
    fitdata = load(strcat(path_fitdata, filenames{i}(1:end-4), "_goright.mat"));
    
    ratnames{i} = fitdata.ratname;
    initpt.loglike(i) = fitdata.loglik;
    initpt.fit_params(i) = sum(fitdata.fit);
    initpt.ntrials(i) = length(fitdata.initial_pt);
    initpt.bic(i) = sum(fitdata.fit)*log(length(fitdata.initial_pt)) - 2*fitdata.loglik;

end


% mean((initpt.bic - lapsemod.bic)./initpt.ntrials)
cols = cbrewer('qual', 'Dark2', 5);

figure();
hb = bar(lapsemod.bic - initpt.bic);
hb(1).FaceColor =  cols(3,:);

ylabel('BIC (inattention) - BIC (motor error)')
xticklabels(ratnames)
xtickangle(60)
box off;

if saving
    set(gcf(), 'Units', 'inches', 'Position', [3,2,6,5]);
    savethisfig(gcf(), path.fig_save + "fig4_rt_modelcomparison")
end


%%
function[ci_bias, ci_lapse, mean_bias, mean_lapse] = compute_lapsebias_CI(prm)


ci_idx = [max(1,floor(length(prm.rc)*0.025)), ceil(length(prm.rc)*0.975)];

% first compute bias CI :
bias_mod = sort(prm.rc(:,4) - prm.lc(:,4));
mean_bias = mean(bias_mod);
ci_bias(1) = mean_bias - bias_mod(ci_idx(1));
ci_bias(2) = bias_mod(ci_idx(2)) - mean_bias;



% now comput lapse CI:
ll = prm.rc(:,1) - prm.lc(:,1);
lr = prm.lc(:,1) + prm.lc(:,2) - (prm.rc(:,1) + prm.rc(:,2));
lapse_mod = sort(ll - lr);
% lapse_mod = sort(prm.lc(:,3)  + prm.lc(:,4) - (prm.rc(:,3) + prm.rc(:,4)));
% lapse_mod = sort(prm.lc(:,4) - prm.rc(:,4) + prm.rc(:,3) - prm.lc(:,3));
mean_lapse = mean(lapse_mod);
ci_lapse(1) = mean_lapse - lapse_mod(ci_idx(1));
ci_lapse(2) = lapse_mod(ci_idx(2)) - mean_lapse;

end





function[] = format_prm_scatter(ax,xmin, xmax)

xlim([xmin, xmax]);
ylim([xmin, xmax]);
plot([xmin, xmax], [xmin, xmax], 'k:')
set(ax, 'linewidth', 2)
set(ax,'TickDir','out')
xlabel('Data');
ylabel('Model');
set(ax, 'FontSize', 14);

axis square;
end