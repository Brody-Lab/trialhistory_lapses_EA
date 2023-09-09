function[boot_prm] = bootstrap_psych_fit(data_master, varargin)

p = inputParser;
addParameter(p,'h_ind', []);
addParameter(p,'pR_ind', []);
addParameter(p,'trial_back', +1);
% +1 means 1 trial after the conditioned trial
addParameter(p, 'xreg','gamma'); % gamma/ Delta
addParameter(p, 'nboots', 1000);
addParameter(p, 'ax', []);
addParameter(p, 'plot_bootstraps', false)
addParameter(p, 'seed', 1)
parse(p,varargin{:});

xreg = p.Results.xreg;
rng(p.Results.seed);

if ~isempty(p.Results.h_ind) & ~isempty(p.Results.pR_ind)
    ids = find(data_master.hits(1:end) == p.Results.h_ind & data_master.pokedR(1:end) == p.Results.pR_ind) + p.Results.trial_back;
    ids(ids<1) = [];
    ids(ids>length(data_master.(xreg))) = [];
elseif ~isempty(p.Results.h_ind)
    ids = find(data_master.hits(1:end) == p.Results.h_ind) + p.Results.trial_back;
    ids(ids<1) = [];
    ids(ids>length(data_master.(xreg))) = [];
elseif isempty(p.Results.h_ind) & isempty(p.Results.pR_ind)
    ids = 1:length(data_master.(xreg));
end

pokedR = data_master.pokedR;
x = data_master.(xreg);
ax = p.Results.ax;

if p.Results.plot_bootstraps
    plotfit = true;
    plotdata = true;
    ploterrorbar = true;
    if isempty(ax)
        figure();
        ax = gca();
        hold on;
    end
else
    plotfit = false;
    plotdata = false;
    ploterrorbar = false;
    ax = [];
    
end

% data_master
boot_prm = zeros(p.Results.nboots, 4);
parfor b = 1:p.Results.nboots  % resample data
    
    id_sample = datasample(ids, length(ids));
    data = [];
    data.pokedR = pokedR(id_sample);
    data.(xreg) = x(id_sample); 
    prm = pbups_psych_gamma(data, ...
        'xreg', p.Results.xreg, ...
        'plotfit', plotfit, ...
        'plotdata', plotdata,...
        'ploterrorbar', ploterrorbar,...
        'fitLineColor', [0.8, 0.8, 0.8],...
        'axHandle', ax);
    drawnow;
    
%     sens = prm.beta(3);
%     bias = prm.beta(4);
%     left_lapse = prm.beta(1);
%     right_lapse = 1 - prm.beta(1) - prm.beta(2);
    
    boot_prm(b,:) = prm.beta;
end


