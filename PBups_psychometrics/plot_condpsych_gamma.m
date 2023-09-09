function[vals, prm] = plot_condpsych_gamma(data, h_ind, pR_ind, alpha, ax, col)

bino_ci_lower = @(n,k) 1-betainv(1-alpha/2, n-k+1, k);
bino_ci_upper = @(n,k) 1-betainv(alpha/2, n-k, k+1);

binc = unique(data.gamma);
id = find(data.hits(1:end-1) == h_ind & data.pokedR(1:end-1) == pR_ind) + 1;

try
    vals(:,2) = grpstats(data.pokedR(id), data.gamma(id));
    n = grpstats(ones(length(id),1), data.gamma(id), 'numel');

    vals(:,1) = vals(:,2) - bino_ci_lower(n,vals(:,2).*n);
    vals(:,3) = bino_ci_upper(n,vals(:,2).*n) - vals(:,2);

    errorbar(ax, binc, vals(:,2), vals(:,1), vals(:,3), 'LineStyle', 'none','LineWidth', 2.5, 'Color', col);
    xlabel('Gamma');
    ylabel('Frac right choices')
    ylim([-0.0 1.]);
catch
    warning('plotting didnt work');
end
[P, ~, exitflag] = fit_logistic4(data.pokedR(id), data.gamma(id));
x_s = linspace(binc(1)-1, binc(end)+1, 100);
plot(x_s, sig4(P,x_s),  'LineStyle', '-', 'LineWidth', 2.5, 'Color', col);

prm.beta = [P.gamma0 P.gamma1 P.sens P.bias];
prm.trust = exitflag>0;
prm.ci = zeros(4,2);
end



function y=sig4(P,x)


y= P.gamma0 + P.gamma1./(1+ exp(-P.sens*(x-P.bias)));


end
