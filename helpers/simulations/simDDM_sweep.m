clear;

%{

This script sweeps through bound and lambda values and saves psychometric paramaters
(conditioned and unconditioned)

Can be easily configured to run free/choice interrogation tasks
Can be easily configured to run for continuous or discrete stimuli
Can be configured to change the drift instead of initial point with history
but the drift setting is not as extensi vely tested.

Check out script analyze_simDDM_param_sweeps to look at analysis code for
plotting and summarizing the simulations

%}

dbstop if error


p = struct();

p.ntrials = 40000;
p.dt = 0.005;

p.seed = 8;
p.saving = true;

p.task_type = "free_choice";
p.stim_type = "continuous";
p.hist_effect = "initial_pt";
p.history_biased = true;

p.mus = -0.5:0.1:0.5;
p.pe = -0.01;
p.te = 0.65;
p.tc = 0.65;
p.pc = 0.2;

p.lambda = -2:0.5:5;
% p.sigma_sens = 0.5:0.5:1.5;
p.sigma_sens = 0.1:0.1:0.6;

p.bound = 0.5:0.1:1.4;


%%
% i = 0;

for bound = p.bound(3:end)
for sigma_sens = p.sigma_sens
    for lambda = p.lambda
        for pc = p.pc
            for tc = p.tc

                i = i+1;
                p.this_pc(i) = pc;
                p.this_tc(i) = tc;
                p.this_bound(i) = bound;
                p.this_sigma_sens(i) = sigma_sens;
                p.this_lambda(i) = lambda;

                fprintf('\n simulating, lambda = %0.2f, sigma_sens = %0.2f, bound = %0.2f, pc = %0.2f, tc = %0.2f', lambda, sigma_sens, bound, pc, tc);

                data = simDDM('bound', bound, 'lambda', lambda, ...
                    "sigma_sens", sigma_sens,...
                    'ntrials', p.ntrials,...
                    'hist_effect', p.hist_effect,...
                    'stim_type', p.stim_type, ...
                    'task_type', p.task_type, ...
                    'history_biased', p.history_biased,...
                    'history_params', make_history_params(pc, -pc, tc, -tc),...
                    'mus', p.mus,...
                    'dt', p.dt,...
                    'seed', p.seed);

                %                     figure()
                clf; hold on;
                ax = gca();
                idx_pc = find(data.hits(1:end-1) == 1)+1;
                data_new.pokedR = data.pokedR(idx_pc);
                data_new.gamma = data.gamma(idx_pc);
                data_new.Delta = data.Delta(idx_pc);
                xreg = 'gamma';
                prm = pbups_psych_gamma(data_new,...
                    'xreg', xreg, ...
                    'plotdata', true,...
                    'ploterrorbar', false,...
                    'fitLineColor', 'k');

                prm_c = plot_conditioned_psychometrics(data, ax, [], 'plot_data', true, 'xreg', xreg, 'binwd_delta', 0.1);
                hline(0, 'k:');
                hline(1, 'k:');
                ylim([-0.05 1.05])
                xlim([-0.8, 0.8]);
                set(ax, 'FontSize', 12);


                xticks([]);
                title(strcat(string(i), ': ', p.task_type, ', ', p.stim_type));
                drawnow;

                sum(abs(data.histbias) > bound)

                %                   figure()
                %                     histogram(data.T)
                p.prm(i) = prm;
                p.prm_c(i) = prm_c;

                %                     fprintf('\n bias mod %d,', p.prm_c(i).beta_rc(4) - p.prm_c(i).beta_lc(4));


            end
        end
    end


end
end




%%
if p.saving == true
filename = strcat('simdata_', p.task_type, '_', p.stim_type, '_bound_lambda_sigmasens_large_pc_med_tc.mat');
save(filename, 'p');
end


%%
% clear;
% load data
master_startup;
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
savethisfig(gcf, path.fig_save + "fig1_bound_lambda_sweeps")
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
savethisfig(gcf, path.fig_save + "fig1_bound_lambda_sweeps_legend")

