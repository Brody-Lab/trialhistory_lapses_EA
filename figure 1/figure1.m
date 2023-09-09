%% generate figure 1
master_startup;


%% plot psych schematics
psych = @(gamma0, gamma1, sens, bias, x) gamma0 + gamma1./(1+exp(-sens*(x-bias)));

x = -55:0.001:55;
sens = 0.15;

%% trial history psychometric schematic

bias = 10;

figure(1001); hold on;
plot(x, (psych(0, 1, sens, bias, x) + psych(0, 1, sens, -bias, x))/2, 'LineWidth', 3, 'Color', c_uncond);
plot(x, psych(0, 1, sens, bias, x), 'LineWidth', 3, 'Color', c_pL);
plot(x, psych(0, 1, sens, -bias, x), 'LineWidth', 3, 'Color', c_pR);
plot([-bias -bias], [-2 2], ':', 'Color', c_pR, 'LineWidth', 2);
plot([bias bias], [-2 2], ':', 'Color',  c_pL, 'LineWidth', 2);
set_plot_psych(gca(), gcf(), path.fig_save + 'fig1_right_biased');


%% lapse schematic

bias = 0;
laps = 0.18;

figure(1002); hold on;
plot(x, (psych(0, 1-laps, sens, bias, x) + psych(laps, 1-laps, sens, -bias, x))/2, 'LineWidth', 3, 'Color', c_uncond);
plot(x, ones(size(x)), ':', 'Color',  'k', 'LineWidth', 2);
set_plot_psych(gca(), gcf(), path.fig_save + 'fig1_lapse_unbiased');


%% plot theoretical psychometrics

uinitpt = 0.45:-0.1:-0.45;
theta = 0.5;
sim_mus = -40:40;
sigma2 = 2;

expf = @(bnd) exp(2*sim_mus*bnd/sigma2);
cols = flip([get_linear_cmap(c_pL, 5); flip(get_linear_cmap(c_pR, 5))]);

figure(1004); clf; hold on;
for i = 1:length(uinitpt)
    
    A = theta - uinitpt(i);
    B = -theta - uinitpt(i);
    Pa(:,i) = (expf(-B) - 1)./(expf(-B) - expf(-A));
    Pa(sim_mus == 0,i) = -B/(A-B);
    plot(sim_mus, Pa(:,i), 'color', cols(i,:), 'LineWidth', 3);
    
end



legend(arrayfun(@num2str, uinitpt, 'UniformOutput', 0),  'location', 'eastoutside');
legend('boxoff');
xlim([-20 20])
set_plot_psych(gca(), gcf(), path.fig_save + 'fig1_initpt_theory', {'Drift rate'}, {'Fraction hit bound B'});





%% plot the PW and PL kernels

pc = 0.34;  % 0.34 0.15
tc = 0.5;   % 0.4 0.82

pe = -0.2;
te = 0.8;

x = [0:20];
offset = 0.3;

figure(1007); clf;
subplot(1,2,1); hold on;
stem(x, pc .* tc.^x, 'filled', 'color', c_pR)
stem(x, -pc .* tc.^x, 'filled', 'color', c_pL)
ylim([-max(abs(pc), abs(pe))- offset, max(abs(pc), abs(pe))+ offset]); 
xlim([-5 25])
yticks([0])

subplot(1,2,2); hold on;
stem(x, -pe .* te.^x, 'filled', 'color', c_pL)
stem(x, pe .* te.^x, 'filled', 'color', c_pR)
ylim([-max(abs(pc), abs(pe))- offset, max(abs(pc), abs(pe))+ offset]); 
xlim([-5 25])
yticks([0])

fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [0.7*fig_pos(3) 0.7*fig_pos(4)];
fig.PaperPositionMode = 'auto';


print('-bestfit', path.fig_save + 'fig1_history_kernels', '-dpdf')





%% plot trajectory of initial points over a whole run 
% N.B. this simulation is not controlled with a seed
% this will be slightly different from the actual figure in the manuscript

bound = 0.6;
mus = -12:0.5:12; %-12
ntrials = 100;

data = simDDM('bound', bound,...
    'ntrials', ntrials,...
    'hist_effect', 'initial_pt',...
    'history_biased', true,...
    'history_params', make_history_params(pc, pe, tc, te),...
    'mus', mus,...
    'seed', 5);


figure(1005); cla;hold on;
cols = [get_linear_cmap(c_pL, 6); flip(get_linear_cmap(c_pR, 6))];

tr_toplot = 20;
bins = discretize(data.histbias, [-0.6, -0.45:0.1:0.45, 0.6]);
ubins = unique(bins(~isnan(bins)));

plot(data.histbias(1:tr_toplot), 'k','LineWidth', 2.5);
axis off
xax = 1:tr_toplot;
for i = 1:length(ubins)
    scatter(xax(bins(1:tr_toplot) == ubins(i)), ...
        data.histbias(bins(1:tr_toplot) == ubins(i)), 350, cols(ubins(i),:),'filled','marker','square')
end
xlim([5.5 tr_toplot+5]);
ylim([-0.8 0.8]);
plot(xlim,[0 0],'k-.')
plot([5.5 5.5],ylim,'k')

% printing trial_type on top
x_offset = 0.3;
y_pos = 0.75;
fontsize = 15;
for i = 6:tr_toplot
    
    if data.pokedR(i-1) == 1
        col = c_pR;
        y = y_pos;
    else
        col = c_pL;
        y = -y_pos;
    end
    
    
    if data.hits(i-1) == 1
        text(i-x_offset, y, 'PW', 'Fontsize', fontsize, 'Color', col)
    else
         text(i-x_offset, y-0.2, 'PL', 'Fontsize', fontsize, 'Color', col)
    end
    
    
end

fig = gcf;
fig_pos = fig.PaperPosition;
fig.PaperSize = [0.7*fig_pos(3) 0.7*fig_pos(4)];
fig.PaperPositionMode = 'auto';

print('-bestfit', path.fig_save + 'fig1_initpt_traj', '-dpdf')



%% plot the apparent lapse psychometrics
figure(1006);clf;

pc = 0.34;  
tc = 0.4;   

pe = -0.01;
te = 0.65;

bound = 0.6;
mus = -12:0.5:12; 
ntrials = 80000;

data = simDDM('bound', bound,...
    'ntrials', ntrials,...
    'hist_effect', 'initial_pt',...
    'history_biased', true,...
    'history_params', make_history_params(pc, pe, tc, te),...
    'mus', mus,...
    'seed', 5);


idx_pc = find(data.hits(1:end-1) == 1)+1;
data_new.pokedR = data.pokedR(idx_pc);
data_new.gamma = data.gamma(idx_pc);
prm = pbups_psych_gamma(data_new,... 
                         'plotdata', false,... 
                         'ploterrorbar', false,... 
                         'fitLineColor', 'k');
hline(1, 'k:');
                 
prm_c = plot_conditioned_psychometrics(data, gca(), [], 'plot_data', false);
xlabel('Drift rate');
ylabel('Fraction chose bound B');
xticks([]);
xticklabels([]);
set_plot_psych(gca(), gcf(), path.fig_save + 'fig1_apparent_lapse', {'Drift rate'}, {'Fraction hit bound B'});