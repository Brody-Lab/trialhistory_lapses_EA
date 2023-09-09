master_startup;

%%
figure(1001); clf;

pe = -0.01;
te = 0.65;

tc_small = 0.01; 
tc_med = 0.65; 
tc_large = 0.84;

pc_small = 0.14;  
pc_large = 0.2;  

pnum = simulate_and_plot(pc_small, pe, tc_small, te, 1);
pnum = simulate_and_plot(pc_large, pe, tc_small, te, pnum);
pnum = simulate_and_plot(pc_small, pe, tc_med, te, pnum);
pnum = simulate_and_plot(pc_large, pe, tc_med, te, pnum);
pnum = simulate_and_plot(pc_small, pe, tc_large, te, pnum);
pnum = simulate_and_plot(pc_large, pe, tc_large, te, pnum);


set(gcf, 'Units', 'inches', 'Position', [3,2,19,10])
savethisfig(gcf(), path.fig_save + "fig1_supp2_timescale_mag_trialhistory")




%% compute functions
function[pnum] = simulate_and_plot(pc, pe, tc, te, i)

bound = 0.725;
mus = -0.8:0.2:0.8;
mus(mus == 0) = [];
ntrials = 40000;

data = simDDM('bound', bound,...
    'ntrials', ntrials,...
    'hist_effect', 'initial_pt',...
    'sigma_sens', 0.3,...
    'dt', 0.005,...
    'history_biased', true,...
    'history_params', make_history_params(pc, pe, tc, te),...
    'mus', mus,...
    'seed', 8);


sum(abs(data.histbias) > bound)


ax = subplot(3,14, [i, i+1,i+2]); hold on;
ntrials_toplot = 50;
padding = 5;
xlim([-padding ntrials_toplot+padding]);
ylim([-bound-0.2 bound+0.2]);
hline(bound, 'k:', '', 1.);
hline(-bound, 'k:','', 1.);
hline(0,'k-')
plot(1:ntrials_toplot , data.histbias(1:ntrials_toplot), 'color',[115,118,140]/255,'LineWidth', 1.5);
text(28, bound, 'bound','BackgroundColor', 'w', 'FontSize', 14) 
text(28, -bound, 'bound','BackgroundColor', 'w', 'FontSize', 14) 
vline(-padding,'k', '', 1.5)
ylabel('Initial point', 'FontSize', 14);
yticks([]);
h = gca; 
h.XAxis.Visible = 'off';
axis square;

ax = subplot(3, 14, [i+4,i+5]); hold on;
idx_pc = find(data.hits(1:end-1) == 1)+1;
data_new.pokedR = data.pokedR(idx_pc);
data_new.gamma = data.gamma(idx_pc);
prm = pbups_psych_gamma(data_new,... 
                         'plotdata', false,... 
                         'ploterrorbar', false,... 
                         'fitLineColor', 'k',...
                         'fitLineWidth', 1);
                     
    
prm_c = plot_conditioned_psychometrics(data, ax, [],...
    'plot_data', false,...
    'fitLineWidth',1.5);
hline(0, 'k:');
hline(1, 'k:');
ylim([-0.05 1.05])
xlim([mus(1), mus(end)]);
set(ax, 'FontSize', 12);

xticks([]);
drawnow;
pnum = i+7;
end

