% wrapper for computing which parameters are significantly modulated by
% trial history:

nboots = 100;

path_pkgdata = '~/ondrive/analysisDG/PBups_trialhistory/data/packaged_reaction_time_data/';
files = dir([path_pkgdata, 'sess_rawdata_*.mat']);
filenames = arrayfun(@(x) x.name,files,'UniformOutput',false);
filenames = filenames(~contains(filenames, 'wholeStim'));


for i = 1:length(filenames)
    load([path_pkgdata, filenames{i}]);
    
    data = [];
    data.pokedR = avgdata.pokedR;
    data.gamma = avgdata.gamma;
    data.hits = avgdata.hits;

    tic()
    boot_rc = bootstrap_psych_fit(data, 'nboots', nboots, 'h_ind', 1, 'pR_ind', 1, 'trial_back', 1);
    boot_lc = bootstrap_psych_fit(data, 'nboots', nboots, 'h_ind', 1, 'pR_ind', 0, 'trial_back', 1);
    toc()
    
    for pn = 1:4
        [p,h] = ranksum(boot_rc(:,pn), boot_lc(:,pn))
    end
    
end

