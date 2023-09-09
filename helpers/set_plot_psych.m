function[] = set_plot_psych(ax, gcf, filename, varargin)


if ~isempty(varargin)
    if length(varargin) == 1
        xlab = varargin{1};
    else
        xlab = 'Stimulus strength';
    end
    if length(varargin) == 2
        xlab = varargin{1};
        ylab = varargin{2};
    end
else
    xlab = 'Stimulus strength';
    ylab = 'Fraction chose right';
end
box off;
yticks([0, 0.5, 1]); ylim([-0.01 1.01]);
xticks([-50, 0, 50]); xticklabels([-100, 0 , 100]);
set(ax, 'linewidth', 2)
set(ax,'TickDir','out')
xlabel(xlab);
ylabel(ylab);
set(ax, 'FontSize', 15);
axis square;

savethisfig(gcf, filename)

end