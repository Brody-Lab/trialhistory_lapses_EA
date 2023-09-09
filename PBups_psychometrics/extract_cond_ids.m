function[ids] = extract_cond_ids(d, cond)
    % this function extracts the ids of the trial for which
    % condition cond is satisfied. session boundaries are respected
    % d is a structure with following fields:
    % d.sessid: vector of session ids corr to each trial 
    % d.hits: 1 hit 0 err 
    % d.pokedR: 1 poked right 0 poked L
 
    
    last_trial_ind = find(diff([d.sessid, 1]));

    switch cond
        case 'last_trial'
            ids = last_trial_ind;
        case 'not_last_trial'
            ids = setdiff(1:length(d.sessid), last_trial_ind);
        case 'corr'
            ids = find(d.hits == 1);
        case 'err'
            ids = find(d.hits == 0);
        case 'post_corr'
            ids = setdiff(find(d.hits == 1), last_trial_ind)+1;
        case 'post_err'
            ids = setdiff(find(d.hits == 0), last_trial_ind)+1;
        case 'post_right_corr'
            ids = setdiff(find(d.hits == 1 & d.pokedR == 1), last_trial_ind)+1;
        case 'post_left_corr'
            ids = setdiff(find(d.hits == 1 & d.pokedR == 0), last_trial_ind)+1;
        case 'post_right_err'
            ids = setdiff(find(d.hits == 0 & d.pokedR == 1), last_trial_ind)+1;
        case 'post_left_err'
            ids = setdiff(find(d.hits == 0 & d.pokedR == 0), last_trial_ind)+1;

           
    end

end
