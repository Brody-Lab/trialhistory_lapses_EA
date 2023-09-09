% check that same number of parameters were fit for each rat

master_startup;


for r = 1:length(ratdata)
    
%     fprintf('===== Processing data for rat %s (%d of %d) ======\n', ratdata(r).name, r, length(ratdata))
    
    for l = 1:length(fits.labels)
        
        label = fits.labels{l};
        load(fullfile(ratdata(r).(label).folder, ratdata(r).(label).filename));
        which_params{l}(r,:) =  fit;
        
        if (l == 3) 
            fprintf('%s, %d \n', ratdata(r).(label).filename, fit(10));
        end
            
        
    end
    
end

for l = 1:length(fits.labels)
    
    subplot(2,2,l)
    imagesc(which_params{l});
    title(fits.labels{l})
    
end