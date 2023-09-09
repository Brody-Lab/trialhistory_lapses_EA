
master_startup;

%%
boot = 0;
nboot = 1;
plotting = 1;

summary = struct();

if plotting
    figure()
end

for r = 1:length(ratdata)
    
    
    fprintf('\n\n===== Processing data for rat %s (%d of %d) ======\n', ratdata(r).name, r, length(ratdata))
    tic()
    
    data = load(fullfile(ratdata(r).datafolder, ratdata(r).dataname), 'avgdata');
    name{r} = ratdata(r).name;
    ntrials(r) = length(data.avgdata.pokedR);
    
    % sanity check
    if ntrials(r) ~= ratdata(r).ntrials
        error("things are not matching up, something is wrong!")
    end
    
    % compute conditioned psychometrics for this base data
    if boot == 1
        rng('shuffle')
        seed = randi(10^8, nboot, 1);
    else
        seed = [];
    end
    
    % compute psychometric parameters of data
    [ll.base(r, :), lr.base(r,:), bias.base(r,:)] = return_condpsych_params(data,...
        boot,...
        nboot,...
        seed,...
        [],...
        plotting);
    if plotting
        title(ratdata(r).name);
        drawnow;
    end
    
    fprintf('\t finished processing base data \n')
    
    % now run for each kind of model
    for l = 1:length(fits.labels)
        
        label = fits.labels{l};
        load(fullfile(ratdata(r).(label).folder, ratdata(r).(label).filename));
        
        % compute psychometric parameters of fits
        [ll.(label)(r, :), lr.(label)(r,:), bias.(label)(r,:)] = return_condpsych_params(data,...
            boot,...
            nboot,...
            seed,...
            p_goright,...
            plotting);
        
        if plotting
            title(strcat(ratdata(r).name, ': ', label ));
            drawnow;
        end
        
        params_ML.(label)(r,:) = ML_params;
        param_names.(label) = name;
        
        % sanity check
        id_ll = label + "_loglik";
        if ratdata(r).(id_ll) ~= loglik
            error("things are not matching up, something is wrong!")
        end
        
        fprintf('\t finished processing model: %s \n', label)
        
    end
    
    
    fprintf('\t okay, finished with this rat took %1f s \n\n', toc())
    
    
    
    
    
end

summary.name = name;
summary.ntrials = ntrials;
summary.ll = ll;
summary.lr = lr;
summary.bias = bias;
summary.params_ML = params_ML;
summary.param_names = param_names;


save(path.mat_save + "summary_most_recent_aug17_2022.mat", 'summary')





%% history of compared fits


% with appropriate prior
% identifier_pbups = {"hist_initpt_jun18", "hist_initpt_lapse_jun18"};
% identifier_bing = {"hist_initpt_jun18", "hist_initpt_lapse_jun18"};

% with appropriate prior and inattention
% identifier_pbups = {"bing_jun18", "hist_initpt_jun18", "hist_initpt_lapse_jun18", "hist_initpt_lapse_inattention_pw_jun18"};
% identifier_bing = {"bing_jun18", "hist_initpt_jun18", "hist_initpt_lapse_jun18", "hist_initpt_lapse_inattention_pw_jun18"};

% no prior and probability matching lapses
% identifier_bing = {"hist_initpt_nov27", "hist_initpt_lapse_nov27"};
% identifier_pbups = {"hist_initpt_jan02", "hist_initpt_lapse_jan10"};

% with slow drift and probability matching lapses?
% identifier_bing = {"hist_initpt_may12", "hist_initpt_lapse_may12"};
% identifier_pbups = {"hist_initpt_may12", "hist_initpt_lapse_may12"};


% with appropriate prior, inattention and initialization of probmatch from
% initial point + fixed implementation of lapse mod_beta
% identifier_pbups = {"bing_jun18", "hist_initpt_jun18", "hist_initpt_lapse_jun18", "hist_initpt_lapse_inattention_pw_jun18"};
% identifier_bing = {"bing_jun18", "hist_initpt_jun18", "hist_initpt_lapse_jun18", "hist_initpt_lapse_inattention_pw_jun18"};
% save('summary_all_jul31.mat')

