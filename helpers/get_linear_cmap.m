function[cols] = get_linear_cmap(c, num_col)

R = linspace(c(1), 0.9, num_col);
G = linspace(c(2), 0.9, num_col);
B = linspace(c(3), 0.9, num_col);

cols = [R' G' B'];
end