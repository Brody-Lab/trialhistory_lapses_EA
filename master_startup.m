clear;

% basepath to this directory
path.base = "/Users/dikshagupta/ondrive/analysisDG/PBups_trialhistory/manuscript/";
     
% relative paths from basepath to the code directory
addpath(genpath(path.base + "figure_code_for_upload/"));

% relative paths from basepath to the figure and sourcedata directory
path.fig_save =  path.base + "figure_pdfs_from_upload/";
path.mat_save =  path.base + "source_data/";
path.data_save =  path.base + "data_for_upload/";


% load processed data with summary statistics
fprintf('Loading ratdata \n')
load(path.mat_save + "valid_ratdata.mat");


%%

% colors for post right and post left psychometric curves
c_postright = get_linear_cmap([0 164 204]./255, 10);
c_postleft = get_linear_cmap([233 115 141]./255, 10);
c_pR = c_postright(1,:);
c_pL = c_postleft(1,:);
c_uncond = [40 51 74]./255;



