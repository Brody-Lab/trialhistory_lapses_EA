function[data] = simDDM(varargin)

%{

This function simulates the DDM trajectorys for given number of trials following Euler-Maruyama method given:  
task type: free_choice or interrogation
stim_type: continuous or discrete
hist_effect: initial_pt, drift or both
bias: i.e. the history bias
bound: value of bound 
bound_crate: collapsing rate of bound
lambda: value of leak
mus: trial difficulties (uniformly sampled)
dt: timestep
seed: for reproducibility
history_biased: whether or not to have history influences
history_params: see make_history_params for format and compute_history_bias

simDDM_one_trial is the major workhorse of this function.

The discrete implementation assumes poisson stimulus. This differs from
Bing's method because there is no adaptation, "bias parameter" or lapse
probability.

If the task_type is interrogation, then trial length is sampled from the
window [0.2 1.0]s with uniform probability.


%}



p = inputParser;
addParameter(p, "task_type", "free_choice", @(x) ismember(x, {'free_choice', 'interrogation'}));
addParameter(p, "stim_type", "continuous", @(x) ismember(x, {'continuous', 'discrete'}));
addParameter(p, "hist_effect", "initial_pt", @(x) ismember(x, {'initial_pt', 'drift', 'both'}));
addParameter(p, "ntrials", 20000);
addParameter(p, "dt", 0.0001);
addParameter(p, "seed", ceil(now));
addParameter(p, "bound", 1);
addParameter(p, "bound_crate", 0.);   % collapsing rate of bound
addParameter(p, "sigma_sens", 1/10);
addParameter(p, "lambda", 0);
addParameter(p, "mus", -8:1:8);
addParameter(p, "history_biased", true);
addParameter(p, "history_params", make_history_params(0.25, -0.25, 0.8, 0.8));
parse(p, varargin{:});

% make new parameters based on inputs
params = fields(p.Results);
for i = 1:length(params)
    eval([params{i} '=p.Results.' params{i} ';' ])
end

rng(seed);

mu = datasample(mus, ntrials);
[choices, outcomes, T, Delta, histbias] = deal(nan(size(mu)));


for tr = 1:ntrials
    
%     if mod(tr, 10000) == 0
%         display(tr)
%     end
    
    if history_biased && tr > 1
        histbias(tr) = compute_history_bias(history_params, choices, outcomes, tr);
    else
        histbias(tr) = 0;
    end
    
    [choices(tr), T(tr), Delta(tr), ~] = simDDM_one_trial(task_type,...
                                               stim_type,...
                                               hist_effect,...
                                               histbias(tr),...
                                               bound,...
                                               bound_crate,...
                                               lambda, ...
                                               sigma_sens, ...
                                               mu(tr),...
                                               mus,...
                                               dt);
%     plot(x);
%     ylim([-bound, bound]);
%     xlim([0 1000]);
%     drawnow;
%     pause(0.5);
    
    
    if mu(tr) == 0
        outcomes(tr) = double(rand()>0.5);
    else
        outcomes(tr) = choices(tr) == (mu(tr) > 0);
    end
        
end


data.pokedR = choices;
data.hits = outcomes;
data.Delta = Delta;
data.gamma = mu;
data.T = T;
data.histbias = histbias;
data.params = p.Results;

end