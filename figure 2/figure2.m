master_startup;


example_rats = {'B052', 'T261', 'B104'};  

%% plotting examples of different types of history effects that are observed in the dataset:


fprintf('\nPlotting example fits \n\n')
figure(2001); clf;

idx = find_rat(example_rats, ratdata);

bin_wd = 3;
range = 32;
markersize = 12;
ylabels = {'Fraction choose right', '', '', ''};

for r = 1:length(idx)
    
    fprintf(' ...processing rat %s\n', ratdata(idx(r)).name)
    
    d = load(fullfile(fullfile(path.data_save, 'pbups'), ratdata(idx(r)).dataname));
    
    h(r) = subplot(1,3,r); hold on;
    
    % compute PC conditioned psychometric curves
    rc = pbups_psych_gamma(d.avgdata, 'xreg', 'Delta',...
        'pR_ind', true,...
        'h_ind', true, ...
        'fitLineColor', c_pR, ...
        'dataLineColor', c_pR, ...
        'errorbarColor', c_pR, ...
        'dataFaceColor', c_pR, ...
        'range_delta', range ,...
        'binwd_delta', bin_wd, ...
        'fitLineWidth', 1., ...
        'dataLineWidth', 0.5,...
        'xlab', '(#R - #L) clicks', ...
        'xticknum', [-20:20:20], ...
        'dataMarkerSize', markersize);
    
    lc = pbups_psych_gamma(d.avgdata, 'xreg', 'Delta',...
        'pR_ind', false,...
        'h_ind', true, ...
        'fitLineColor', c_pL, ...
        'dataLineColor', c_pL, ....
        'errorbarColor', c_pL, ...
        'dataFaceColor', c_pL, ...
        'range_delta', range,...
        'binwd_delta', bin_wd, ...
        'fitLineWidth', 1., ...
        'dataLineWidth', 0.5,...
        'xlab', '(#R - #L) clicks', ...
        'xticknum', [-20:20:20], ...
        'dataMarkerSize', markersize);
    
    
    % compute PC psychometric curve
    idx_pc = find(d.avgdata.hits(1:end-1) == 1) + 1;
    data.pokedR = d.avgdata.pokedR(idx_pc);
    data.Delta = d.avgdata.Delta(idx_pc);
    base = pbups_psych_gamma(data, 'xreg', 'Delta', ...
        'range_delta', range,...
        'binwd_delta', bin_wd, ...
        'fitLineWidth', 0.8, ...
        'dataLineWidth', 0.5,...
        'dataMarkerSize', markersize,...
        'ylab', ylabels{r}, ...
        'xlab', '(#R - #L) clicks', ...
        'xticknum', [-20:20:20], ...
        'axisFontSize', 15);
    
    h(r).XTickLabelRotation = 0;
    title(sprintf('Example rat %d', r))
    
    
    
end

set(gcf, 'Units', 'inches', 'Position', [3,2,15,5])
savethisfig(gcf(), path.fig_save + "fig2_examplerats")



%% plot population summary across rats


d = reorganize_prm(ratdata, 'prm');
params = {'sens', 'bias', 'left_lapse', 'right_lapse'};
titles = {'Sensitivity', 'Threshold', 'Left lapse rate', 'Right lapse rate'};
fprintf('\n ==== \n')


idx = find_rat(example_rats, ratdata);

for i = 1:length(params)
    figure(2002); clf; hold on;

%     subplot(2,2,i); hold on;
    
    plotdata = [d.rc.(params{i})' d.lc.(params{i})'];
    plot_param_scatters(plotdata, gca(), c_pR, c_pL);
    ylabel(titles{i})
    
    % highlight example rats
    plot(gca(), [1 1.5], [d.rc.(params{i})(idx)' d.lc.(params{i})(idx)']', '-o',...
        'markerfacecolor', [0.8 0.8 0.8],...
        'markersize',6,...
        'LineWidth', 1, ...
        'color',  [0.4 0.4 0.4]);
    
    % plotting mean again, so it is on top
    plot([1 1.5], median(plotdata), '-o',...
        'markerfacecolor', 'k',...
        'markersize',6,...
        'LineWidth', 1.5, ...
        'color',  'k');
    scatter(1, median(plotdata(:,1)), 50, c_pR, 'filled');
    scatter(1.5, median(plotdata(:,2)), 50, c_pL, 'filled');
    
    set(gca(), 'linewidth', 1)
    
    [p,h] = signrank(d.lc.(params{i})', d.rc.(params{i})');
    
    if h == 0
        fprintf('\n %s is NOT significantly different with p-value %d\n', titles{i}, p)
    else
        fprintf('\n %s IS significantly different with p-value %d\n', titles{i}, p)
    end
    
    if i == 1
        ylim([0 0.5]);
    elseif i == 2
        ylim([-8, 8]);
    elseif i == 3 | i == 4
        ylim([0 0.35]);
    end
    
    
    set(gcf, 'Units', 'inches', 'Position', [3,2,3,4])
    fig = gcf();
    savethisfig(gcf(), path.fig_save + "fig2_paramscatters" + string(i))

    
end







