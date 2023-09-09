
function[ratdata] = get_valid_fits(ratdata, fits)


for r = 1:length(ratdata)
    
    fprintf('Processing rat %s (%d of %d)\n', ratdata(r).name, r, length(ratdata))
    
    
    ratdata(r).valid = true;
    
    data  = load(fullfile(ratdata(r).datafolder,ratdata(r).dataname), 'rawdata');
    ntrials = length(data.rawdata);
    ratdata(r).ntrials = ntrials;
    switch ratdata(r).dataset
        case 'bing'
            sigma = nan(ntrials, 1);
            for tr = 1:ntrials
                sigma(tr) = length(data.rawdata(tr).leftbups) + length(data.rawdata(tr).rightbups);
            end
        case 'pbups'
            sigma = [data.rawdata.Sigma];
    end
    
    if sum(sigma == 0)/ntrials > fits.min_totalclickszero
        ratdata(r).valid = false;
    end
    
    if ntrials < fits.min_ntrials
        ratdata(r).valid = false;
    end
    
    
    for l = 1:length(fits.labels)
        
        identifier = fits.identifier.(ratdata(r).dataset){l};
        id_ll = fits.labels{l} + "_loglik";
        id_ML = fits.labels{l} + "_MLparams";

        switch ratdata(r).dataset
            case 'bing'
                files = dir(strcat(fits.fit_dir, 'bing_fits/chrono_', ratdata(r).name, '*', identifier,'.mat'));
            case 'pbups'
                files = dir(strcat(fits.fit_dir, 'pbups_fits/chrono_rawdata_',ratdata(r).name,'*',identifier,'.mat'));
        end
        
        if ~isempty(files)
            
            % get the file with the best loglikelihood
            ll_temp = nan(length(files),1);
            for f = 1:length(files)
                load(strcat(files(f).folder, filesep, files(f).name));
                ntrials = length(p_goright);
                ll_temp(f) = loglik;
                if exp(ll_temp(f)/ntrials) > 1
                    ll_temp(f) = nan;
                end
            end
            [logll, file_num] = nanmax(ll_temp);
            
            % make sure the likelihood is valid, sometimes because of bad
            % initialization the fit fails
            if isnan(exp(logll/ntrials))
                warning(identifier + " :only bad fit exists");
                ratdata(r).(fits.labels{l}).filename = nan;
                ratdata(r).(fits.labels{l}).folder = nan;
                ratdata(r).(id_ll) = nan;
                ratdata(r).valid = false;
                ratdata(r).(id_ML) = nan;
                
            else
                
                % assign the info
                ratdata(r).(fits.labels{l}).filename = files(file_num).name;
                ratdata(r).(fits.labels{l}).folder = files(file_num).folder;
                ratdata(r).(id_ll) = logll;
                ratdata(r).(id_ML) = ML_params;
                if (exp(logll/ntrials) < fits.min_ll) | (exp(logll/ntrials) > 0.98)
                    ratdata(r).valid = false;
                end
                
            end
            
        else
            warning(identifier + " :fit file doesn't exist");
            % otherwise mark that fit doesn't exist
            ratdata(r).(fits.labels{l}).filename = nan;
            ratdata(r).(fits.labels{l}).folder = nan;
            ratdata(r).(id_ll) = nan;
            ratdata(r).valid = false;
            ratdata(r).(id_ML) = nan;
        end
        
    end
    
end


% only return valid rats
ratdata = ratdata([ratdata(:).valid] == 1);


end


