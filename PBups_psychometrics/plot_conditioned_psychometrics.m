function[prm] = plot_conditioned_psychometrics(data, ax_postcorr, ax_posterr, varargin)
% data.pokedR, data.gamma, data.hits
% axes for post corr, post err
% trans is unfortunately, whether or not to plot in transparency true/false

p = inputParser;

addParameter(p,'xreg', 'gamma');
addParameter(p,'plot_data', []);
addParameter(p, 'compute_fit', true);
addParameter(p, 'plot_fit', true);
addParameter(p, 'trans', false);
addParameter(p, 'fitLineWidth', []);
addParameter(p, 'fit_pokedR', []);
addParameter(p, 'binwd_delta', 2);
addParameter(p, 'axisFontSize', 18);
addParameter(p, 'plot_fitpokedR', false); % plot the model fits without fitting them

parse(p,varargin{:});

c_postright = get_linear_cmap([0 164 204]./255, 10);
c_postleft = get_linear_cmap([233 115 141]./255, 10);

if p.Results.trans
    c_pR = c_postright(8,:);
    c_pL = c_postleft(8,:);
    lw = 1.;
    pt = true;
else
    c_pR = c_postright(1,:);
    c_pL = c_postleft(1,:);
    lw = 2;
    pt = true;
    
end
[c_pR_fit, c_pR_data] = deal(c_pR);
[c_pL_fit, c_pL_data] = deal(c_pL);

if ~isempty(p.Results.fit_pokedR)
    pt = true;
    c_pR_data = c_postright(1,:);
    c_pL_data = c_postleft(1,:);
    c_pR_fit = c_postright(1,:);
    c_pL_fit = c_postleft(1,:);
    lw = 2;
end

if ~isempty(p.Results.plot_data)
    pt = p.Results.plot_data;
end

if ~isempty(p.Results.fitLineWidth)
    lw = p.Results.fitLineWidth;
end

% POST CORRECT psychometrics
prm_rc = pbups_psych_gamma(data,  'xreg', p.Results.xreg, ...
    'fit_pokedR', p.Results.fit_pokedR, ...
    'pR_ind', true, 'h_ind', true, ...
    'compute_fit', p.Results.compute_fit,...
    'plotfit', p.Results.plot_fit,...
    'plotdata', pt, 'ploterrorbar', pt,...
    'axHandle', ax_postcorr, ...
    'fitLineWidth', lw, ...
    'axisFontSize', p.Results.axisFontSize, ...
    'binwd_delta', p.Results.binwd_delta, ...
    'fitLineColor', c_pR_fit, ...
    'dataLineColor', c_pR_data, ...
    'errorbarColor', c_pR_data, ...
    'dataFaceColor', c_pR_data, ...
    'plot_fitpokedR', p.Results.plot_fitpokedR);

prm_lc = pbups_psych_gamma(data, 'xreg', p.Results.xreg, ...
    'fit_pokedR', p.Results.fit_pokedR, ...
    'pR_ind', false, 'h_ind', true, ...
    'compute_fit', p.Results.compute_fit,...
    'plotfit', p.Results.plot_fit,...
    'plotdata', pt, 'ploterrorbar', pt,...
    'axHandle', ax_postcorr, ...
    'fitLineWidth', lw, ...
    'axisFontSize', p.Results.axisFontSize, ...
    'binwd_delta', p.Results.binwd_delta, ...
    'fitLineColor', c_pL_fit, ...
    'dataLineColor', c_pL_data, ....
    'errorbarColor', c_pL_data, ...
    'dataFaceColor', c_pL_data, ...
    'plot_fitpokedR', p.Results.plot_fitpokedR);

if ~isempty(ax_posterr)
    % POST ERROR psychometrics
    prm_re = pbups_psych_gamma(data, 'xreg', p.Results.xreg, ...
        'fit_pokedR', p.Results.fit_pokedR, ...
        'plot_fitpokedR', p.Results.plot_fitpokedR,...
        'pR_ind', true, 'h_ind', false, ...
        'compute_fit', p.Results.compute_fit,...
        'plot_fit', p.Results.plot_fit,...
        'plotdata', pt, 'ploterrorbar', pt,...
        'axHandle', ax_posterr, ...
        'fitLineWidth', lw, ...
        'axisFontSize', p.Results.axisFontSize, ...
        'binwd_delta', p.Results.binwd_delta, ...
        'fitLineColor', c_pR_fit, ...
        'dataLineColor', c_pR_data, ...
        'errorbarColor', c_pR_data, ...
        'dataFaceColor', c_pR_data);
    
    prm_le = pbups_psych_gamma(data, 'xreg', p.Results.xreg, ...
        'fit_pokedR', p.Results.fit_pokedR, ...
        'plot_fitpokedR', p.Results.plot_fitpokedR,...
        'pR_ind', false, 'h_ind', false, ...
        'compute_fit', p.Results.compute_fit,...
        'plot_fit', p.Results.plot_fit,...
        'plotdata', pt, 'ploterrorbar', pt,...
        'axHandle', ax_posterr, ...
        'binwd_delta', p.Results.binwd_delta, ...
        'fitLineWidth', lw, ...
        'axisFontSize', p.Results.axisFontSize, ...
        'fitLineColor', c_pL_fit, ...
        'dataLineColor', c_pL_data, ...
        'errorbarColor', c_pL_data, ...
        'dataFaceColor', c_pL_data);
end


if ~isempty(prm_rc) & p.Results.plot_fitpokedR
    prm.fracR_rc = prm_rc.fracR;
    prm.fracR_lc = prm_lc.fracR;
    prm.fitfracR_rc = prm_rc.fitfracR;
    prm.fitfracR_lc = prm_lc.fitfracR;
    prm.ntrials_rc = prm_rc.ntrials;
    prm.ntrials_lc = prm_lc.ntrials;
    if ~isempty(ax_posterr)
        prm.fracR_re = prm_re.fracR;
        prm.fracR_le = prm_le.fracR;
        prm.fitfracR_re = prm_re.fitfracR;
        prm.fitfracR_le = prm_le.fitfracR;
        prm.ntrials_re = prm_re.ntrials;
        prm.ntrials_le = prm_le.ntrials;
    end
end

if ~isempty(prm_rc)
    prm.beta_rc = prm_rc.beta;
    prm.beta_lc = prm_lc.beta;
    prm.trust = prm_rc.trust & prm_lc.trust;
    if ~isempty(ax_posterr)
        prm.beta_re = prm_re.beta;
        prm.beta_le = prm_le.beta;
    else
        [prm.beta_re,prm.beta_le] = deal([0 0 0 0]);
    end
else
    [prm.beta_rc,prm.beta_lc,prm.beta_re,prm.beta_le] = deal([0 0 0 0]);
end


end