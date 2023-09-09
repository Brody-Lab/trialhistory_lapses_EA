
function [merge_s] = append_struct(aa_s,bb_s)

%APPENDSTRUCT horizontally appends two structures 

% Convert structures to tables
aa_t = struct2table(aa_s);
bb_t = struct2table(bb_s);

% Concatonate tables
merge_t = [aa_t; bb_t];


% Convert table to structure
merge_s = table2struct(merge_t);

end
