function[choice, T, Delta, x] = simDDM_one_trial(task_type, stim_type, hist_effect, bias, bound, bound_crate, lambda, sigma_sens, mu, mus, dt)

%{

This function simulates the DDM trajectory for 1 trial following Euler-Maruyama method given:
task type: free_choice or interrogation
stim_type: continuous or discrete
hist_effect: initial_pt, drift or both
bias: i.e. the history bias
bound: value of bound
lambda: value of leak
mu: trial difficulty
dt: timestep

simDDM is the major workhorse of this function

The discrete implementation assumes poisson stimulus. This differs from
Bing's method because there is no adaptation, "bias parameter" or lapse
probability.

If the task_type is interrogation, then trial length is sampled from the
window [0.2 1.0]s with uniform probability.

%}

click_R = 1;
click_L = -1;

switch task_type
    case "free_choice"
        T = 20;
        x = nan(T/dt,1);
    case "interrogation"
        T = round(0.2+ 0.8*rand(),3);
        x = nan(round(T/dt), 1);
end

switch hist_effect
    case "initial_pt"
        x(1) = bias;
        
    case "drift"
        x(1) = 0;
        switch stim_type
            case "continous"
                mu = mu + bias;
            case "discrete"
                click_R = click_R + bias;
                click_L = click_L + bias;
        end
        
    case "both"
        error("Not implemented properly yet")
        x(1) = bias;
        mu = mu + bias;
end


switch stim_type
    case "continuous"
        
        samples = 0;
        Delta = 0;
        
        for i = 2:length(x)
            
            %             samples = samples + mu*dt + sqrt(dt)*randn();
            
            sample = mu*dt + sigma_sens*sqrt(dt)*randn();
            
            % add N(mu dt, sig^2 dt) to x((i-1) dt)
%             if lambda == 0
                
%                 pu = zeros(length(mus),1);
%                 for m = 1:length(mus)
%                     pu(m) = 2/length(mus)*exp(-(mus(m)*dt - sample)^2/(2*dt));
%                 end
                pu = 2/length(mus) .* exp(-(mus.*dt - sample).^2./(2*sigma_sens*dt));
                
                num_right = log(sum(pu(mus>0)));
                num_left = log(sum(pu(mus<0)));
                llr = num_right - num_left;
%                 x(i) = (x(i-1)*exp(dt*lambda)) + llr;
                x(i) = x(i-1) + x(i-1)*lambda*dt + llr;

                %                 x(i) = 0;
                %                 norm = 0;
                %                 for m = 1:length(mus)
                %                     a_m = 1/length(mus)*exp(samples*mus(m) - 0.5*(i-1)*dt*mus(m)^2);
                %                     b_m = 1/length(mus)*exp(-samples*mus(m) - 0.5*(i-1)*dt*mus(m)^2);
                %                     x(i) = x(i) + a_m;
                %                     norm = norm + a_m + b_m;
                %                 end
                %                 x(i) = x(i)/norm;
                %                 x(i) = log(x(i)/(1- x(i))) + x(1);
%             else
%                 warning('Not implemented yet');
%                 x(i) = mu * dt + (x(i-1) * exp(dt*lambda)) + (1*sqrt(dt) *randn());
%             end
            B = 2* bound./(1+exp(bound_crate*(i-1)));
            if abs(x(i)) >= B
                x(i) = B * sign(x(i));
                break;
            end
        end
        
    case "discrete"
        % first generate discrete stimulus
        Rbar = (40*exp(mu))./(1+exp(mu));
        R = cumsum(exprnd(1/Rbar,40,1));
        L = cumsum(exprnd(1/(40-Rbar),40,1));
        rightbups = R(R <= T);
        leftbups = L(L <= T);
        Delta = length(rightbups) - length(leftbups);
        
        % and then run the simulation
        for i = 2:length(x)
            
            % decay past value
            if lambda == 0
                x(i) = x(i-1);
            else
                x(i) = x(i-1)*exp(dt*lambda);
            end
            
            % add leftbups
            while(~isempty(leftbups) && leftbups(1) <= i*dt)
                x(i) = x(i) + click_L + randn()*sigma_sens;
                leftbups = leftbups(2:end);
            end
            
            % add rightbups
            while(~isempty(rightbups) && rightbups(1) <= i*dt)
                x(i) = x(i) + click_R + randn()*sigma_sens;
                rightbups = rightbups(2:end);
            end
            
            B = 2* bound./(1+exp(bound_crate*(i-1)));
            if abs(x(i)) >= B
                x(i) = B * sign(x(i));
                break;
            end
        end
        
end

choice = x(i) > 0;
if strcmp(task_type, "free_choice")
    T = i*dt;
end

end