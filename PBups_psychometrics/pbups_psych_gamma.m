function[prm] = pbups_psych_gamma(data, varargin)

p = inputParser;

addParameter(p,'h_ind', []);
addParameter(p,'pR_ind', []);
addParameter(p,'trial_back', +1);
% +1 means 1 trial after the conditioned trial

addParameter(p, 'xreg','gamma'); % gamma/ Delta
addParameter(p, 'binwd_delta', 2); % binwidth for deltaclicks discretization - (if normalizing by sqrt total clicks, 0.5)
addParameter(p, 'range_delta', 40); % assumed to be symmetric - (if normalizing by sqrt total clicks, 5)

addParameter(p, 'fit_pokedR', []); % fit predicted responses with a pscyhometric curve, requires 0/1
addParameter(p, 'fit_gp', false); % fit a GP model to interpolate. require p(R)
addParameter(p, 'plot_fitpokedR', false); % plot the model fits without fitting them

addParameter(p,'compute_fit',true);
addParameter(p,'plotfit', true);
addParameter(p,'plotdata', true);
addParameter(p,'ploterrorbar', true);
addParameter(p,'axHandle', []);
addParameter(p,'grid', 'off');

addParameter(p,'errorbarColor',[0 0 0]);
addParameter(p,'dataLineColor', [0 0 0]);
addParameter(p,'dataFaceColor', [0 0 0]);
addParameter(p,'dataMarkerSize', 18);
addParameter(p,'dataLineStyle', '.');
addParameter(p,'dataLineWidth', 1);
addParameter(p,'fitLineStyle', '-');
addParameter(p,'fitLineWidth', 2.5);
addParameter(p,'fitLineColor',[0 0 0]);

addParameter(p,'axisFontSize', 18);
addParameter(p,'xtickangle', 60);
addParameter(p,'xlab', []);
addParameter(p,'ylab','Fraction chose right');
addParameter(p,'xticknum',-4.5:1.5:4.5);
addParameter(p,'yticknum',[0, 0.5, 1]);

addParameter(p, 'bootstrap_seed', []);   % if this is nonempty a bootstrap sample is drawn with its value as the seed

parse(p,varargin{:});


plot_fitpokedR = p.Results.plot_fitpokedR;
xregname = p.Results.xreg;

xlab_dict.gamma = 'Stimulus strength';
xlab_dict.Delta = 'Delta clicks';
xticknum_dict.gamma = [-4.5:1.5:4.5];
xticknum_dict.Delta = [-20:20:20];
% if normalizing by sqrt clicks: use th
% xticknum_dict.Delta = [-8:2:8];



if p.Results.plotfit || p.Results.plotdata || p.Results.ploterrorbar
    if ~isempty(p.Results.axHandle)
        ax = p.Results.axHandle;
        hold on;
    else
        ax = gca();
        hold on;
    end
    
    if isempty(p.Results.xlab)
        xlab = xlab_dict.(xregname);
    else
        xlab = p.Results.xlab;
    end
    
    if isempty(p.Results.xlab)
        xticknum = xticknum_dict.(xregname);
    else
        xticknum = p.Results.xticknum;
    end
end

if ~isempty(p.Results.h_ind) & ~isempty(p.Results.pR_ind)
    ids = find(data.hits(1:end) == p.Results.h_ind & data.pokedR(1:end) == p.Results.pR_ind) + p.Results.trial_back;
    ids(ids<1) = [];
    ids(ids>length(data.(xregname))) = [];
else
    ids = 1:length(data.(xregname));
    
end

xreg = data.(xregname)(ids);
pokedR = data.pokedR(ids);

% sampling for bootstrapping
if ~isempty(p.Results.bootstrap_seed)
    rng(p.Results.bootstrap_seed)
    [xreg, boo] = datasample(xreg, length(xreg));
    pokedR = pokedR(boo);
end

if ~isempty(p.Results.fit_pokedR)
    fit_pokedR = p.Results.fit_pokedR(ids)';
    if ~isempty(p.Results.bootstrap_seed)
        fit_pokedR = fit_pokedR(boo);
    end
end



if  p.Results.compute_fit
    if isempty(p.Results.fit_pokedR)
        [P, ~, exitflag] = fit_logistic4(pokedR, xreg);
        prm.beta = [P.gamma0 P.gamma1 P.sens P.bias];
        prm.trust = exitflag>0;
        prm.ci = zeros(4,2);
    elseif ~isempty(p.Results.fit_pokedR) & p.Results.fit_gp
        gp = fitrgp(xreg', fit_pokedR);
        prm = [];
    elseif ~isempty(p.Results.fit_pokedR)
        [P, ~, exitflag] = fit_logistic4(fit_pokedR, xreg);
        prm.beta = [P.gamma0 P.gamma1 P.sens P.bias];
        prm.trust = exitflag>0;
        prm.ci = zeros(4,2);
    end
    
    if p.Results.plotfit
        switch xregname
            case 'gamma'
                x_s = linspace(min(xreg)-1, max(xreg)+1, 100);
            case 'Delta'
                x_s = linspace(-p.Results.range_delta-8, p.Results.range_delta+8, 100);
        end
        if ~isempty(p.Results.fit_pokedR) & p.Results.fit_gp
            
            plot(ax, x_s, predict(gp, x_s'),  p.Results.fitLineStyle,...
                'LineWidth', p.Results.fitLineWidth, ...
                'color', p.Results.fitLineColor);
        elseif ~plot_fitpokedR
            plot(ax, x_s, sig4(P,x_s),  p.Results.fitLineStyle,...
                'LineWidth', p.Results.fitLineWidth, ...
                'color', p.Results.fitLineColor);
        end
    end
    
else
    prm = [];
    
end

alpha = 0.05;
bino_ci_lower = @(n,k) 1-betainv(1-alpha/2, n-k+1, k+eps());  % eps for robustness to nans when n=k or k=0
bino_ci_upper = @(n,k) 1-betainv(alpha/2, n-k+eps(), k+1);
% https://sigmazone.com/binomial-confidence-intervals/

switch xregname
    case 'gamma'
        uxreg = unique(xreg);
        for i = 1:length(uxreg)
            ntrials(i) = sum(xreg == uxreg(i));
            fracR(i) = sum(pokedR(xreg == uxreg(i)))./ntrials(i);
        end
    case 'Delta'
        edges = -p.Results.range_delta - 0.25 : p.Results.binwd_delta : p.Results.range_delta + 0.25;
        uxreg = round(edges(1:end-1) +  p.Results.binwd_delta/2);
        %%
        %         if normalizing by sqrt total clicks use this
        %                 edges = -p.Results.range_delta - 0.05 : p.Results.binwd_delta : p.Results.range_delta + 0.05;
        %                 uxreg = edges(1:end-1) +  p.Results.binwd_delta/2;
        %%
        for i = 1:length(uxreg)
            ntrials(i) = sum(xreg < edges(i+1) & xreg > edges(i));
            fracR(i) = sum(pokedR(xreg < edges(i+1) & xreg > edges(i)))./ntrials(i);
            if ~isempty(p.Results.fit_pokedR) & plot_fitpokedR
                fit_fracR(i) = sum(fit_pokedR(xreg < edges(i+1) & xreg > edges(i)))./ntrials(i);
            end
        end
end

if p.Results.plotdata || p.Results.ploterrorbar
    %     dotsize = (((ntrials-min(ntrials))/max(ntrials))*15).^1.5 + 15;
    if p.Results.ploterrorbar
        cifracRneg = fracR - bino_ci_lower(ntrials, fracR.*ntrials);
        cifracRpos = bino_ci_upper(ntrials,fracR.*ntrials) - fracR;
        errorbar(ax, uxreg, fracR, cifracRneg, cifracRpos, ...
            'LineStyle', 'none',...
            'color',p.Results.errorbarColor,...
            'markerfacecolor',p.Results.dataFaceColor,...
            'linewidth',p.Results.dataLineWidth,...
            'markersize',p.Results.dataMarkerSize,'capsize',0);
    end
    plot(ax, uxreg, fracR, p.Results.dataLineStyle,...
        'markerfacecolor',p.Results.dataFaceColor,...
        'markersize',p.Results.dataMarkerSize,...
        'LineWidth', p.Results.dataLineWidth, ...
        'color', p.Results.dataLineColor);
end

if plot_fitpokedR
    cifracRneg = fit_fracR - bino_ci_lower(ntrials, fit_fracR.*ntrials);
    cifracRpos = bino_ci_upper(ntrials,fit_fracR.*ntrials) - fit_fracR;
    ebregionplot(uxreg, fit_fracR, cifracRneg, cifracRpos, ...
        p.Results.fitLineColor, 'ax', ax);
    prm.fracR = fracR;
    prm.fitfracR = fit_fracR;
    prm.ntrials = ntrials;
end



if p.Results.plotfit || p.Results.plotdata || p.Results.ploterrorbar
    grid(p.Results.grid)
    xlabel(ax, xlab); ylabel(ax,p.Results.ylab);
    xticks(ax, xticknum); xtickangle(ax,p.Results.xtickangle);
    yticks(ax,p.Results.yticknum);
    ylim(ax,[-0.01 1.01]);
    set(ax, 'linewidth', 1.)
    set(ax,'TickDir','out')
    set(ax, 'FontSize', p.Results.axisFontSize);
    box off;
    axis(ax, 'square');
end

end

function y=sig4(P,x)


y= P.gamma0 + P.gamma1./(1+ exp(-P.sens*(x-P.bias)));


end
