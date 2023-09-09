master_startup;
nrats = length(ratdata);


%% Win-stay lose-switch across rats - ignoring session boundaries, effects would be negligible

p_repeat_win = nan(nrats,1);
p_repeat_loss = nan(nrats,1);
ci_win = nan(nrats,2);
ci_loss = nan(nrats,2);

for r = 1:nrats
    
    fprintf('Processing rat %s\n', ratdata(r).name);
    load(fullfile(fullfile(path.data_save, 'pbups'), ratdata(r).dataname));
    
    idx_win = find(avgdata.hits(1:end-1) == 1);
    num_repeats_win = sum(avgdata.pokedR(idx_win) == avgdata.pokedR(idx_win+1));
    p_repeat_win(r) = num_repeats_win/length(idx_win);
    [ci_win(r,:),~,~] = bino_confidence(length(idx_win), num_repeats_win, 0.95);
    
    idx_loss = find(avgdata.hits(1:end-1) == 0);
    num_repeats_loss = sum(avgdata.pokedR(idx_loss) == avgdata.pokedR(idx_loss+1));
    p_repeat_loss(r) = num_repeats_loss/length(idx_loss);
    [ci_loss(r,:),~,~] = bino_confidence(length(idx_loss), num_repeats_loss, 0.95);
    
end

%%
figure(2006); clf; hold on;

ax = gca();
cols = cbrewer('qual', 'Set1', nrats);

for r = 1:nrats
    errorbar(p_repeat_win(r), p_repeat_loss(r),...
        ci_loss(r,1) - p_repeat_loss(r),...
        ci_loss(r,2) - p_repeat_loss(r),...
        ci_win(r,1) - p_repeat_win(r),...
        ci_win(r,2) - p_repeat_win(r),...
        'd',...
        'MarkerFaceColor', cols(r,:),...
        'MarkerEdgeColor', cols(r,:),...
        'Color', cols(r,:),...
        'MarkerSize', 2,...
        'capsize',0,...
        'linewidth',0.2)
end

ylim([0.3,0.58])
xlim([0.4,0.68])
hline(0.5, 'k:')
vline(0.5,'k:')
ylabel('P(repeat) post loss');
xlabel('P(repeat) post win');
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;
axis(ax, 'square');

set(gcf, 'Units', 'inches', 'Position', [3,2,5,5])


savethisfig(gcf(), path.fig_save + "fig2_supp_winstayloseswitch")


%% plot bias-lapse modulation scatters


[ci_bias, ci_lapse] = deal(nan(length(ratdata),2));
[mean_bias, mean_lapse, mean_uncond_lapse] = deal(nan(length(ratdata),1));


for i = 1:length(ratdata)
    
    display(ratdata(i).name);
    
    prm = ratdata(i).boot_prm;
    
    % compute 95% CIs
    ci_idx = [max(1,floor(length(prm.rc)*0.025)), ceil(length(prm.rc)*0.975)];
    
    % first compute bias CI :
    bias_mod = sort(prm.rc(:,4) - prm.lc(:,4));
    mean_bias(i) = mean(bias_mod);
    ci_bias(i,1) = mean_bias(i) - bias_mod(ci_idx(1));
    ci_bias(i,2) = bias_mod(ci_idx(2)) - mean_bias(i);
    
    % now comput lapse CI:
    ll = prm.rc(:,1) - prm.lc(:,1);
    lr = prm.lc(:,1) + prm.lc(:,2) - (prm.rc(:,1) + prm.rc(:,2));
    lapse_mod = sort(ll - lr);
    mean_lapse(i) = mean(lapse_mod);
    ci_lapse(i,1) = mean_lapse(i) - lapse_mod(ci_idx(1));
    ci_lapse(i,2) = lapse_mod(ci_idx(2)) - mean_lapse(i);
    
    % compute total lapse modulation in the psychometric
    mean_uncond_lapse(i) = mean(1 - prm.uncond(:,2));
    
    
end

bias_0 = (- ci_bias(:,1) + mean_bias) < 0 & (ci_bias(:,2)+ mean_bias) > 0;
lapse_0 = (- ci_lapse(:,1) + mean_lapse) < 0 & (ci_lapse(:,2)+ mean_lapse) > 0;

%%
figure(2003); clf; hold on;
ax = gca();
cols = cbrewer('qual', 'Set3', 5);
plot_error_scatter(ax, mean_bias, mean_lapse, ci_bias, ci_lapse, 1:length(bias_0), cols(4,:),2);

% highlight example rats
example_rats = {'B052', 'T261', 'B104'};   %K214
idx = find_rat(example_rats, ratdata);
plot_error_scatter(ax, mean_bias, mean_lapse, ci_bias, ci_lapse, idx, 'k', 4);
xlim([-12, 6]);
ylim([-0.1, 0.35]);

hline(0, 'k:')
vline(0,'k:')
ylabel('Lapse rate modulation');
xlabel('Threshold modulation');
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 10);
box off;
axis(ax, 'square');

set(gcf, 'Units', 'inches', 'Position', [3,2,5,5])
savethisfig(gcf(), path.fig_save + "fig2_supp_biasmod_scatters")




%% plot a histogram of fraction of lapse that is history modualted

figure(2005); clf; hold on;

subplot(2,2,4); hold on;
ax = gca();
histogram(abs((mean_lapse./2))./mean_uncond_lapse,...
    'binwidth', 0.05,...
    'LineWidth', 1.,...
    'FaceColor', [0.6, 0.6, 0.6]);
xlim([-0.05,0.8]);

xlabel({'Frac of lapse rates','modulated by history'});
ylabel('Number of rats');
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;



% Plot number of trials and stability of accuracy


% Distribution of number of trials
num_bins = 10;

subplot(2,2,1); hold on;
ax = gca();
histogram([ratdata(:).ntrials],...
    'NumBins', num_bins,...
    'LineWidth', 1.,...
    'FaceColor', [0.6, 0.6, 0.6]);
% xlim([-0.05,0.8]);
xticks([10000, 50000, 100000, 150000]);
xlabel('Number of trials');
ylabel('Number of rats');
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;

%%

% Plot stability of accuracy over this period
filt_size = 200;
accuracy_slopes = nan(nrats,1);
accuracy_intercept = nan(nrats, 1);

for r = 1:nrats
    
    fprintf('Processing rat %s\n', ratdata(r).name);
    data = load(fullfile(fullfile(path.data_save, 'pbups'), ratdata(r).dataname));
    smooth_data = smooth(data.avgdata.hits, filt_size);
    
    % do a linear fit
    p = polyfit(1:length(smooth_data),smooth_data,1);
    accuracy_slopes(r) = p(1);
    accuracy_intercept(r) = p(2);
    
end
%%
% accuracy slope
subplot(2,2,2); hold on;
ax = gca();
histogram(accuracy_slopes*10000,...
    'LineWidth', 1.,...
    'NumBins', num_bins,...
    'FaceColor', [0.6, 0.6, 0.6]);
xlabel({'Average change in accuracy','every 10000 trials'});
set(ax, 'linewidth', 1.)
ylabel('Number of rats');
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;

% accuracy intercept
subplot(2,2,3); hold on;
ax = gca();
histogram(accuracy_intercept,...
    'LineWidth', 1.,...
    'NumBins', num_bins,...
    'FaceColor', [0.6, 0.6, 0.6]);
xlabel('Mean accuracy');
ylabel('Number of rats');

set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;


%% save!
set(gcf, 'Units', 'inches', 'Position', [3,2,12,6])
savethisfig(gcf(), path.fig_save + "fig2_supp_rat_ntrials_accuracy_statistics")


%% bias-lapse modulation
function plot_error_scatter(ax, mean_bias, mean_lapse, ci_bias, ci_lapse, ind, col, markersize)

errorbar(mean_bias(ind), mean_lapse(ind),...
    ci_lapse(ind,1),...
    ci_lapse(ind,2),...
    ci_bias(ind,1),...
    ci_bias(ind,2),...
    'd',...
    'MarkerSize', markersize,...
    'capsize',0,...
    'MarkerFaceColor', col,...
    'MarkerEdgeColor', col,...
    'Color', col,...
    'linewidth',0.2)

end