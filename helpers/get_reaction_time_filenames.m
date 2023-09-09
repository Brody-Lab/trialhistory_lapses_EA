function[filenames] = get_reaction_time_filenames(path_pkgdata)


files = dir(path_pkgdata + 'sess_rawdata_*.mat');
filenames = arrayfun(@(x) x.name,files,'UniformOutput',false);
filenames = filenames(~contains(filenames, 'wholeStim'));

end