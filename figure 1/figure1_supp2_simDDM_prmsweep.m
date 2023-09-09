master_startup;

%%

load(path.mat_save + "simdata_free_choice_continuous_bound_lambda_sigmasens_large_pc_med_tc.mat");

% extract post right correct parameters from the struct
for i = 1:length(p.prm)
    lc_left_lapse(i) = p.prm_c(i).beta_lc(1);
    lc_right_lapse(i) = p.prm_c(i).beta_lc(1) + p.prm_c(i).beta_lc(2);
    lc_bias(i) = p.prm_c(i).beta_lc(4);
    rc_left_lapse(i) = p.prm_c(i).beta_rc(1);
    rc_right_lapse(i) = p.prm_c(i).beta_rc(1) + p.prm_c(i).beta_rc(2);
    rc_bias(i) = p.prm_c(i).beta_rc(4);
    
    sens(i) = p.prm_c(i).beta_rc(3);
end

left_lapse = rc_left_lapse - lc_left_lapse;
right_lapse = rc_right_lapse - lc_right_lapse;
bias =  rc_bias - lc_bias;

unique_sigma_sens = unique(p.sigma_sens);
unique_bound = unique(p.this_bound);
unique_bound = unique_bound(1:2:end-1);
unique_lambda = unique(p.this_lambda);

cols = cbrewer('div', 'Spectral', length(unique_sigma_sens));
size = 2:2:200;

figure(1010); hold on; clf;

for b = 1:length(unique_bound)
    idx_bound = p.this_bound == unique_bound(b);
    subplot(1,3, b); hold on;
    
    k = 1;
    for lambda = unique_lambda(1:end-1)
        idx_lambda = p.this_lambda == lambda;
        
        for i = 1:length(unique_sigma_sens)
            
            idx = p.this_sigma_sens == unique_sigma_sens(i);
            this_idx = idx & idx_bound & idx_lambda;
            scatter(bias(this_idx),...
                (right_lapse(this_idx) + left_lapse(this_idx))/2,...
                size(k),...
                'MarkerFaceColor', cols(i,:), ...
                'MarkerEdgeColor', cols(i,:));
        end
        k = k+1;
    end
    xlim([-0.55, 0.05]);
    ylim([-0.01, 0.6]);
    hline(0,'k--');
    vline(0, 'k--');
    axis square;
    set(gca, 'TickDir', 'out', 'TickLength', [0.02 ,0.25]);
    set(gca, 'linewidth', 0.9)
    set(gca, 'FontSize', 13);
    if b > 1
        yticks([]);
    end
    
end

%%
saving = 1;
if saving
    set(0,'DefaultFigureWindowStyle', 'normal')
    set(gcf, 'Units', 'inches', 'Position', [3,2,8,3])
    savethisfig(gcf, path.fig_save + "fig1_supp_bound_lambda_sweeps")
end


%% creat a fake legend plot

figure();
hold on;

k = 1;
for l = 1:length(unique_lambda(1:end-1))
    for i = 1:length(unique_sigma_sens)
        i
        scatter(1,l, size(k), 'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor', 'k');
    end
    k = k+1;
end

for i = 1:length(unique_sigma_sens)
    scatter(1.2,i, 'MarkerFaceColor', cols(i,:), ...
        'MarkerEdgeColor', cols(i,:));
end

xlim([0.8,1.4])
set(gcf, 'Units', 'inches', 'Position', [3,2,8,3])
savethisfig(gcf, path.fig_save + "fig1_supp_bound_lambda_sweeps_legend")

