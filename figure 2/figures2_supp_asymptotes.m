% this script first generates fake psychometric data with lapse rate, number of trials, range and discretization of diifficulties.
% then it fits a psychometric curve and reports the lapse rate
% it also estimates the lapse rate based on the end point at the highest trial difficulty
% and then reports how good the recovery is
master_startup;

figure(); hold on;


default.true_lapse = 0.2;
default.sensitivity = 0.125;
default.bias = 0;
default.ntrials = 10000;
default.stimrange = 40;
default.stimspacing = 1;
default.lapse_range = 0:0.05:0.35;

ntrials = 2000:2000:20000;
stimrange = 10:5:45;


%%
% plot ntrials, bias = 0
toplot = compute_data('ntrials', ntrials, default);
plot_scatters(toplot, ntrials, 'fitlapse', 1, 'Greens');
plot_scatters(toplot, ntrials, 'asymptote', 2, 'Greens');
make_legend(ntrials, 9, 'Greens')

% plot stimrange, bias = 0
toplot = compute_data('stimrange', stimrange, default);
plot_scatters(toplot, stimrange, 'fitlapse', 5, 'Purples');
plot_scatters(toplot, stimrange, 'asymptote', 6, 'Purples');
make_legend(stimrange, 10, 'Purples')

% make psychometric curves for each of the cases - can also be used for
% demonstrating how things are being measured
default.bias = 0;
data =  simulate_data(default,'true_lapse', 0.2);
gca = subplot(3,4,11); cla()
data_temp.pokedR = data.choice;
data_temp.Delta = data.stimulus;
prm = pbups_psych_gamma(data_temp,...
    'xreg', 'Delta', ...
    'plotfit', true,...
    'plotdata', true,...
    'binwd_delta', 4,...
    'axisFontSize', 14,...
    'ploterrorbar', true);
        

% plot ntrials, bias = 10
default.bias = 10;
toplot = compute_data('ntrials', ntrials, default);
plot_scatters(toplot, ntrials, 'fitlapse', 3, 'Greens');
plot_scatters(toplot, ntrials, 'asymptote', 4, 'Greens');

% plot stimrange, bias = 10
toplot = compute_data('stimrange', stimrange, default);
plot_scatters(toplot, stimrange, 'fitlapse', 7, 'Purples');
plot_scatters(toplot, stimrange, 'asymptote', 8, 'Purples');

data =  simulate_data(default,'true_lapse', 0.2);
gca = subplot(3,4,12); cla();
data_temp.pokedR = data.choice;
data_temp.Delta = data.stimulus;
prm = pbups_psych_gamma(data_temp,...
    'xreg', 'Delta', ...
    'plotfit', true,...
    'binwd_delta', 4,...
    'plotdata', true,...
    'axisFontSize', 14,...
    'ploterrorbar', true);

%% save the figure 
set(gcf, 'Units', 'inches', 'Position', [3,2,18,12])
savethisfig(gcf(), path.fig_save + "fig2_supp_asymptotes_recovery")
        
        
%%

function[toplot] = compute_data(this_varname, this_variable, default)

toplot = struct();
lapse_range = default.lapse_range;
toplot.lapse_range = lapse_range;

total = length(this_variable) * length(lapse_range);

for vr = 1:length(this_variable)
    for lp = 1:length(toplot.lapse_range)
        
        fprintf('%d remaining \n', total);
        data =  simulate_data(default,...
            'true_lapse', toplot.lapse_range(lp),...
            this_varname, this_variable(vr));
        
        toplot.(this_varname)(vr,lp) = this_variable(vr);
        toplot.truelapse(vr,lp) = lapse_range(lp);
        toplot.fitlapse(vr,lp) = fit_lapse(data);
        toplot.asymptote(vr,lp) = compute_asymptote(data);

        total = total - 1;
    end
end

end



function[] = plot_scatters(toplot, variable, whichlapse, subplot_num, cmap)

s = 45;
colormap = cbrewer('seq', cmap, length(variable));

subplot(3,4,subplot_num); hold on;
for i = 1:length(variable)
    scatter(toplot.truelapse(i,:),...
        toplot.(whichlapse)(i,:),...
        s*ones(size(toplot.truelapse(i,:))),...
        'filled',...
        'MarkerFaceColor', colormap(i,:));
end

% plot(toplot.lapse_range, toplot.lapse_range, 'k');
axis square;

ax = gca();
xmin = toplot.lapse_range(1) - 0.1;
xmax = toplot.lapse_range(end) + 0.1;
xlim([xmin, xmax]);
ylim([xmin, xmax]);
plot([xmin, xmax], [xmin, xmax], 'k:')
set(ax, 'linewidth', 1.5)
set(ax,'TickDir','out')
% xlabel('True lapse rate');
% ylabel('Inferred');
set(ax, 'FontSize', 14);


end



function[] = make_legend(variable, subplot_num, cmap)

subplot(3,4,subplot_num); hold on;
colormap = cbrewer('seq', cmap, length(variable));
for i = 1:numel(variable)
   bubleg(i) = scatter(0,0,'filled', 'MarkerFaceColor',colormap(i,:)); 
   set(bubleg(i),'visible','off')
   legentry{i} = num2str(variable(i));
end
legend(legentry)
legend boxoff
h = gca;
h.XAxis.Visible = 'off';
h.YAxis.Visible = 'off';

end



function[fitlapse] = fit_lapse(data)

data_temp.pokedR = data.choice;
data_temp.Delta = data.stimulus;
prm = pbups_psych_gamma(data_temp,...
    'xreg', 'Delta', ...
    'plotfit', false,...
    'plotdata', false,...
    'ploterrorbar', false);
fitlapse = 1 - prm.beta(2);

end


function[asymptote] = compute_asymptote(data)

stimrange = [max(data.stimulus), min(data.stimulus)];
for sr = 1:length(stimrange)
    idx = data.stimulus == stimrange(sr);
    rchoice(sr) = mean(data.choice(idx));
    if sr == 1
        rchoice(sr) = 1-rchoice(sr);
    end
end
asymptote = sum(rchoice);

end


function[data] = simulate_data(default, varargin)

p = inputParser;
addParameter(p, 'true_lapse', default.true_lapse);
addParameter(p, 'sensitivity', default.sensitivity);
addParameter(p, 'bias', default.bias);
addParameter(p, 'ntrials', default.ntrials);
addParameter(p, 'stimrange', default.stimrange);
addParameter(p, 'stimspacing', default.stimspacing);

parse(p, varargin{:});

data = struct();
data.true_lapse = p.Results.true_lapse;
data.sensitivity = p.Results.sensitivity;
data.bias = p.Results.bias;
data.ntrials = p.Results.ntrials;
data.stimrange = p.Results.stimrange;
data.stimspacing = p.Results.stimspacing;

stims = -data.stimrange:data.stimspacing:data.stimrange+mod(data.stimrange, data.stimspacing);
data.stimulus = randsample(stims, data.ntrials,1);
for i = 1:data.ntrials
    if rand() < data.true_lapse
        data.choice(i) = datasample([0,1], 1);
    else
        p = 1/(1+exp(-data.sensitivity*(data.stimulus(i)-data.bias)));
        data.choice(i) = rand() < p;
    end
end



end
