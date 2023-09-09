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


p = struct();

p.ntrials = 40000;
p.dt = 0.005;

p.seed = 8;
p.saving = false;

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
p.sigma_sens = 0.1:0.1:0.6;

p.bound = 0.5:0.1:1.4;


%%


for bound = p.bound
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
                    p.prm(i) = prm;
                    p.prm_c(i) = prm_c;
                    
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

