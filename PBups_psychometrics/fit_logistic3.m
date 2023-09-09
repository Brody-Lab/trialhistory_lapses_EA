function [P varargout] = fit_logistic3(pokedR, xregressor)


if size(pokedR) ~= size(xregressor)
    error("one choice for each regressor, something is wrong!")
end
x_uniq = unique(xregressor);

nPokedR = [];
nPokes = [];

for c = 1:numel(x_uniq)
    
    idx = x_uniq(c) == xregressor;
    
    nPokedR(c,1) = sum(pokedR(idx)); 
    nPokes(c,1) = sum(idx);
end

% *** bounds ***
lowerbound = [0, 0, -range(xregressor)]; % lapse [0,1], sensitivity [0 5], bias [-40, 40]
upperbound = [1, 20,  range(xregressor)];
% lowerbound = [];
% upperbound = [];


% *** initial parameters ***
param0 = [0.1, 0.5, 0]; % no lapse, mediocre sensitivity, no bias

warning('off', 'MATLAB:nchoosek:LargeCoefficient')

% ***************** %
%       Fmincon     %
% ***************** %
nLL0 = negLogLike_logistic3(param0, nPokes, nPokedR, x_uniq);
if isinf(nLL0) || isnan(nLL0)
    error('Failure to evaluate initial values')
end

[x_fmin, f_fmin, exitflag, output, ~, grad, hessian] = ...
    fmincon(@(param) negLogLike_logistic3(param, nPokes, nPokedR, x_uniq), ...
            param0, ...
            [], [], [], [], ...
            lowerbound, upperbound, ...
            [], ...
            optimset('DiffMinChange', 0.00001, ...
                     'MaxIter', 1000, ...
                     'Algorithm', 'interior-point', ...
                     'Display', 'off')     );
P = struct;
P.lapse  = x_fmin(1);
P.sens   = x_fmin(2);
P.bias   = x_fmin(3);

% *** vargout ***
varargout{1} = f_fmin;
varargout{2} = exitflag;
varargout{3} = output;
varargout{4} = grad;
varargout{5} = hessian;
end