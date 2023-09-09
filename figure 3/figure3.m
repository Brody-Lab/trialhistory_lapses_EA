master_startup;

load(path.mat_save + "summary_most_recent_aug17_2022.mat")
saving = 1;

% compute number of params
nparams.nohistory = 10;
nparams.initpt = 14;
nparams.inattention = 14;

fits.labels = {"nohistory", "initpt", "inattention"};

% lets first augment summary with logllpertrial, bic, aic etc.
for l = 1:length(fits.labels)
    idx_ll = fits.labels{l} + '_loglik';
    label_ll = [ratdata(:).(idx_ll)];
    label_ntrials = [ratdata(:).ntrials];
    label_nparams = nparams.(fits.labels{l});
    summary.logllpertrial.(fits.labels{l}) = exp(label_ll./label_ntrials);
    summary.aic.(fits.labels{l}) = 2*label_nparams - 2*label_ll;
    summary.bic.(fits.labels{l}) = label_nparams*log(label_ntrials) - 2*label_ll;
end

% adding a field for total lapse variation
flds = fieldnames(summary.bias);
for l = 1:length(flds)
    summary.lapse.(flds{l}) = (summary.ll.(flds{l}) - summary.lr.(flds{l}));
end


num_rats = length(ratdata);



%% plotting figure for model comparison (hist vs no hist)


figure('Name', 'model comparion'); hold on;

subplot(2,1,1);
idx = find_winning_model(summary, 1);
bar([histcounts(idx)],...
    'LineWidth', 1.,...
    'FaceColor', [0.6, 0.6, 0.6]);
% title('Winning model (individuals)');
xticks([1,2]);
xticklabels({'No hist + Motor err', 'Hist + Motor err'})
ylabel('Number of rats')
xtickangle(35)
ax = gca();
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 14.5);
box off;

subplot(2,1,2); hold on;
S = summary;
BIC = [S.bic.nohistory'./S.ntrials', S.bic.initpt'./S.ntrials'];
meanBIC = sum(BIC)/length(ratdata);
sem = std(BIC)/sqrt(length(ratdata));
bar(meanBIC,...
    'LineWidth', 1.,...
    'FaceColor', [0.6, 0.6, 0.6]);

er = errorbar([1,2],meanBIC,sem);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

% title('Winning model (population)');
% ylim([min(meanBIC-sem)-100, max(meanBIC+sem) + 10])
ylim([min(meanBIC-sem)-0.005, max(meanBIC+sem) + 0.005])
xticks([1,2]);
xticklabels({'No hist + Motor err', 'Hist + Motor err'})
ylabel('Mean per trial BIC per rat')
xtickangle(35)
ax = gca();
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 14.5);
box off;

if saving == 1
    set(gcf, 'Units', 'inches', 'Position', [3,2,3.5,10])
    savethisfig(gcf(), path.fig_save + "fig3_modelcomparison_hist_nohist")
end

[p,h] = signrank(BIC(:,1), BIC(:,2),'tail','right')
[h,p] = ttest(BIC(:,1), BIC(:,2),'tail','right')



%% plotting figure for model comparison (all models)


figure('Name', 'model comparion'); hold on;

subplot(2,1,1);
idx = find_winning_model(summary,2);
bar([histcounts(idx)],...
    'LineWidth', 1.,...
    'FaceColor', [0.6, 0.6, 0.6]);
% title('Winning model (individuals)');
xticks([1,2,3]);
xticklabels({'No hist', 'Motor err', 'Inattention'})
ylabel('Number of rats')
xtickangle(35)
ax = gca();
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;

subplot(2,1,2);
S = summary;
BIC = [S.bic.nohistory', S.bic.initpt', S.bic.inattention'];
meanBIC = sum(BIC)/length(ratdata);
bar(meanBIC,...
    'LineWidth', 1.,...
    'FaceColor', [0.6, 0.6, 0.6]);
% title('Winning model (population)');
ylim([min(meanBIC)-100, max(meanBIC) + 10])
xticks([1,2,3]);
xticklabels({'No hist', 'Motor err', 'Inattention'})
ylabel('Mean BIC per rat')
xtickangle(35)
ax = gca();
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;

if saving == 1
    set(gcf, 'Units', 'inches', 'Position', [3,2,4,8])
    savethisfig(gcf(), path.fig_save + "fig3_supp_modelcomparison")
end

%% 
figure()
imagesc(BIC' - min(BIC'))
caxis([0, 200])
colorbar()
yticks([1,2,3]);
yticklabels({'No hist', 'Motor error', 'Inattention'})
title('\Delta BIC from the lowest score')
xlabel('Rats')
colormap('copper')
ax = gca();
set(ax, 'linewidth', 1.)
set(ax,'TickDir','out')
set(ax, 'FontSize', 12);
box off;
if saving == 1
    set(gcf, 'Units', 'inches', 'Position', [3,2,8,3])
    savethisfig(gcf(), path.fig_save + "fig3_supp_modelcomparison_perrat")
end

%%

figure('Name', 'Param scatters')
models = {'initpt'};
idx = find_winning_model(summary,1);
cols = cbrewer('seq','Blues',8);
variables = {'bias', 'lapse'};
for var  = 1:length(variables)
    currvar = variables{var};
    [val, val_sort] = get_sort_for_variable(summary, currvar);
    
    for i = 1:length(models)
        subplot(length(models),length(variables), var + (i-1)*length(variables)); hold on;
        plot_from_summary(summary, currvar, models(i), gca())
        ylabel('Model');
    end
    xlabel('Data');
    title(currvar);
end


if saving == 1
    set(gcf, 'Units', 'inches', 'Position', [3,2,8,8])
    savethisfig(gcf(), path.fig_save + "fig3_param_scatters")
end


%% compute R2

% To get R2 and slope values for lapses
mdl = fitlm(S.lapse.base, S.lapse.initpt, 'RobustOpts', 'on', 'intercept', false);
fprintf("\n\nInitpt - Lapse vs base: \n Rsq = %s \n Slope = %s",...
    string(mdl.Rsquared.Adjusted),...
    string(mdl.Coefficients.Estimate));

% To get R2 and slope values for biases
mdl = fitlm(S.bias.base, S.bias.initpt, 'RobustOpts', 'on', 'intercept', false);
fprintf("\n\nInitpt - Bias vs base: \n Rsq = %s \n Slope = %s",...
    string(mdl.Rsquared.Adjusted),...
    string(mdl.Coefficients.Estimate));


%%
function[idx] = find_winning_model(S, compare_type)

if compare_type == 1
    BIC = [S.bic.nohistory', S.bic.initpt'];
    [val, idx] = min(BIC, [], 2);  
elseif compare_type == 2
    BIC = [S.bic.nohistory', S.bic.initpt', S.bic.inattention'];
    [val, idx] = min(BIC, [], 2);
end

end



function[val, val_sort] = get_sort_for_variable(S, varb)
[val, val_sort] = sort(S.(varb).base);
if strcmp(varb, 'bias')
    val_sort = flip(val_sort);
    val = flip(val);
end
end



function[] = format_barplot(b)
cols = cbrewer('seq','Blues',8);

for j = 1:4
    for i = 1:3
        b(i).CData(j,:) = cols(i+3,:);
    end
end
box off;
set(gca, 'FontSize', 20)
set(gca, 'linewidth', 1.5)
set(gca,'TickDir','out')
xticks([1:4]);
yticks([0:25:100])
ylim([0 85])
set(gca,  'xticklabels',{'No hist', 'Motor err', 'Inatt max', 'Inatt match'})
xtickangle(30);

end



function[] = plot_from_summary(S, variable, labels, ax)

a = .4;
dotsize = 80*([S.ntrials]/max([S.ntrials]));

fldnames = fields(S.(variable));
for l = fldnames(2:end)
    xmin = min(vertcat(S.(variable).base, S.(variable).(l{1})), [], 'all');
    xmax = max(vertcat(S.(variable).base, S.(variable).(l{1})), [], 'all');
end
x = xmax - 0.3*(xmax-xmin);
y = xmin + 0.1*(xmax-xmin);

combn = nchoosek(cellstr(horzcat(["base"], labels)), 2);
if isempty(ax)
    figure()
    if length(combn) < 3
        k = 1;
    else
        k = ceil(length(combn)/3);
    end
end

% make the dotsize legend
legend_spacing = 8;
bubsizes = unique(round(dotsize/legend_spacing)*legend_spacing)';
for ind = 1:numel(bubsizes)
   bubleg(ind) = plot(0,0,'o','markersize', sqrt(bubsizes(ind)), 'color',[0 0.4470 0.7410] ); %,'filled',...
%        'MarkerFaceColor', [0 0.4470 0.7410],...
%        'MarkerFaceAlpha', a,...
%        'MarkerEdgeAlpha', a);
   set(bubleg(ind),'visible','off')
   legentry{ind} = num2str(bubsizes(ind)/80*max(S.ntrials));
end
legend(legentry)
legend boxoff

for c = 1:size(combn,1)
    if isempty(ax)
        subplot(k, min(length(combn), 3), c); hold on;
    end
    scatter(S.(variable).(combn{c,1}), S.(variable).(combn{c,2}), dotsize, 'filled',...
        'MarkerFaceColor', [0 0.4470 0.7410],...
        'MarkerFaceAlpha', a,...
        'MarkerEdgeAlpha', a)
    
    format_prm_scatter(gca(), xmin, xmax+0.01)
    xlabel(combn{c,1})
    ylabel(combn{c,2})
    model_combo_label = strcat(combn{c,1}, '_vs_', combn{c,2});
%         text(x,y,sprintf('R: %.2f', S.R2.(model_combo_label).(variable)),'FontSize',16);
end




end


function[] = format_prm_scatter(ax,xmin, xmax)

xlim([xmin, xmax]);
ylim([xmin, xmax]);
plot([xmin, xmax], [xmin, xmax], 'k:')
set(ax, 'linewidth', 1.5)
set(ax,'TickDir','out')
xlabel('Data');
ylabel('Model fit');
set(ax, 'FontSize', 14);

axis square;
end



