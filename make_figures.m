

warning('Edit paths in master_startup before getting started!');

%% generate figures
master_startup;
set(0,'DefaultFigureWindowStyle','normal')


fprintf('\n\n\nGENERATING FIGURE 1 PANELS\n\n\n');
figure1;

fprintf('\n\n\nGENERATING FIGURE 2 PANELS\n\n\n');
figure2;

fprintf('\n\n\nGENERATING FIGURE 3 PANELS\n\n\n');
figure3;

fprintf('\n\n\nGENERATING FIGURE 4 PANELS\n\n\n');
figure4_rtdata_metarat;
figure4_rthistorybias_theory_and_data;
figure4_plot_example_rat_fits;
figure4_modelresults;
