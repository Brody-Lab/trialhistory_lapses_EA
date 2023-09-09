% population summary across rats
function[] = plot_param_scatters(data, ax, c_pR, c_pL)

plot(ax, [1 1.5], data', '-o',...
    'markerfacecolor', [0.8 0.8 0.8],...
    'markersize',2,...
    'LineWidth', 0.5, ...
    'color',  [0.8 0.8 0.8]);

plot(ax, [1 1.5], median(data), '-o',...
    'markerfacecolor', 'k',...
    'markersize',8,...
    'LineWidth', 2.5, ...
    'color',  'k');

scatter(ax, 1, median(data(:,1)), 50, c_pR, 'filled');
scatter(ax, 1.5, median(data(:,2)), 50, c_pL, 'filled');

set(ax, 'linewidth', 1)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;
xlim(ax, [0.75 1.75])
% ylim(ax_se, [0.8 1.5])
xticks(ax, [1 1.5]);
xticklabels(ax, {'','',});
xtickangle(ax, 0)

end
