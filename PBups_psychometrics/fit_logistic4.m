function [P varargout] = fit_logistic4(pokedR, xregressor)

rng shuffle
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
lowerbound = [0, 0, 0.0001*range(xregressor), -range(xregressor)]; % lapse [0,1], sensitivity [0 5], bias [-40, 40]
upperbound = [1.0, 1.0, 100*range(xregressor), range(xregressor)];


function [c,ceq] = boundcon(x)
    c(1) = x(1) + x(2) - 1; % gamm0 + gamm1  <= 1
    c(2) = -(x(1) + x(2));  % gamm0 + gamm1 >= 0
    ceq = [];
end

% *** initial parameters ***
param0 = [0.01, 0.9, 0.01, 0]; % low lapse, mediocre sensitivity, no bias

warning('off', 'MATLAB:nchoosek:LargeCoefficient')

% ***************** %
%       Fmincon     %
% ***************** %
nLL0 = negLogLike_logistic4(param0, nPokes, nPokedR, x_uniq);
if isinf(nLL0) || isnan(nLL0)
    error('Failure to evaluate initial values')
end

[x_fmin, f_fmin, exitflag, output, ~, grad, hessian] = ...
    fmincon(@(param) negLogLike_logistic4(param, nPokes, nPokedR, x_uniq), ...
            param0, ...
            [], [], [], [], ...
            lowerbound, upperbound, ...
            @boundcon, ...
            optimset('DiffMinChange', 1e-12, ...
                     'MaxIter', 5000, ...
                     'Algorithm', 'interior-point', ...
                     'Display', 'off')     );
P = struct;
P.gamma0 = x_fmin(1);
P.gamma1 = x_fmin(2);
P.sens   = x_fmin(3);
P.bias   = x_fmin(4);

% *** vargout ***
varargout{1} = f_fmin;
varargout{2} = exitflag;
varargout{3} = output;
varargout{4} = grad;
varargout{5} = hessian;
end