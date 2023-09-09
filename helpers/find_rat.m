function[idx] = find_rat(ratlist, ratdata)

% find the index of a rat from a ratlist in ratdata

idx = nan(length(ratlist),1);
for i = 1:length(ratlist)
    idx(i) = find(strcmp({ratdata(:).name},ratlist{i}));
end

end
