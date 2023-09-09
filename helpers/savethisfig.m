function[] = savethisfig(fig, filename)
% fig: gcf


fig_pos = fig.PaperPosition;
fig.PaperSize = [0.7*fig_pos(3) 0.7*fig_pos(4)];
fig.PaperPositionMode = 'auto';

print('-bestfit', filename, '-dpdf')

fprintf('YESS QUEEN saved! loc: %s\n', pwd())

end