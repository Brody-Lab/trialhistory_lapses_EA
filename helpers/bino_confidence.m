function [ci, lower, upper] = bino_confidence(n, k, varargin)
input_parser = inputParser;
addOptional(input_parser, 'confidence_level', 0.95, @(x) isnumeric(x) && isscalar(x))
parse(input_parser, varargin{:});
for param = fields(input_parser.Results)'
    eval([param{:} ' = input_parser.Results.' param{:} ';']);
end

% 2016_01_09
% http://www.sigmazone.com/binomial_confidence_interval.htm
%
%   n                   number of samples (successes + failures)
%   k                   number of successes
%   confidence_level    confidence level; takes a value within (0, 1).

n = n(:);
k = k(:);

assert(numel(n) == numel(k), '"n" and "k" need to have the same size')
num = numel(n);

ci = nan(num,2);

% lower confidence bound
ci(:,1) = 1 - betainv( (1+confidence_level)/2, ...
                        n - k + 1, ...
                        k );

% upper confidence bound
ci(:,2) = 1 - betainv( (1-confidence_level)/2, ...
                        n - k, ...
                        k + 1 );

lower = k./n - ci(:,1);
upper = ci(:,2) - k./n;

