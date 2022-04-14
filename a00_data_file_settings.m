%%% do not include clear, clc, or close all command here

% this file provides paths and tree information that is required by
% multiple scripts. It is loaded elsewhere to remove redundant coding. Be
% sure to change the various paths and tree_names to match your
% application.

%% define paths

%%% specify the working path (i.e., where this script is stored)
if ispc==1
    work_path = 'C:\Users\Mark Raleigh\Google Drive\Research\Colorado\Sway_SnowIntercept\MATLAB\';
else
    work_path = '/home/mark/Research/Colorado/Sway_SnowIntercept/MATLAB/';
end

%%% specify the paths where the raw data are stored for each sensor. Note
%%% the order here must match the order in tree_names (below).
%%% NOTE: this code assumes each of these folders has sub-folders each
%%% specific to a collection date with folder name YYYYMMDD. If this is not
%%% the case, you will need to adapt the code a bit.
path_raw_data = {'/home/mark/Research/Colorado/Niwot/Canopy Interception Study/GCDC_Data/01_Main_Tree_Trunk/', ...
             '/home/mark/Research/Colorado/Niwot/Canopy Interception Study/GCDC_Data/03_Small_Tree_Trunk/', ...
             '/home/mark/Research/SnowEx/Tree_Accel/Colorado/SBB/North_Tree/', ...
             '/home/mark/Research/SnowEx/Tree_Accel/Colorado/SBB/South_Tree/', ...
             '/home/mark/Research/SnowEx/Tree_Accel/Colorado/Grand_Mesa/GM08/', ...
             '/home/mark/Research/SnowEx/Tree_Accel/Colorado/Grand_Mesa/GM12A/', ...
             '/home/mark/Research/SnowEx/Tree_Accel/Colorado/Grand_Mesa/GM12B/'};

         
%%% specify the path where you want to store the outputs for the compiled
%%% raw datasets in netCDF format:
if ispc==1
    path_out_comp = 'C:\Users\Mark Raleigh\Big_Data\Colorado\Sway_SnowIntercept\';
else
    path_out_comp = '/home/mark/Research/Colorado/Sway_SnowIntercept/MATLAB/Data/';
end

%%% specify the path where you want to store the outputs for the frequency
%%% analysis (.mat files):
if ispc==1
    path_out_freq = 'C:\Users\Mark Raleigh\Google Drive\Research\Colorado\Sway_SnowIntercept\MATLAB\Data\';
else
    path_out_freq = '/home/mark/Research/Colorado/Sway_SnowIntercept/MATLAB/Data/';
end

%% define tree/sensor characteristics

%%% specify unique names for each site/tree/sensor. The number of entries
%%% must match the number of entries in data_path (above)
tree_names = {'Niwot_Tree01', 'Niwot_Tree02', 'SBB_North', 'SBB_South', 'GM08', 'GM12A', 'GM12B'};

%% tree/sensor specific parameters

%%% maximum and minimum sway frequency (limits)
%%%  (each column corresponds to the tree_name in
%%% same position)
fsway_max = [1.38 0.91];
fsway_min = [0.65 0.38];


%% filtering and smoothing parameters

%%% power threshold (only retain points with power threshold below this)
pthresh = 0.05;

%%% smoothing parameter for splines
SmoothParam_sway = 0.99;     % 0.9 works decently for sway frequency
SmoothParam_Temp = 0.9999;  % close to 1 works best for Temp 

%%% max/min temp
temp_max = 30;
temp_min = -30;

%% visualization

%%% max length of gap (hrs) to show spline
gap_max = 24;

%% checks

%%% check to make sure all sensors have a path and a name
if numel(path_raw_data)~=numel(tree_names)
    error('all sensors must include both a data_path entry and a tree_name entry')
end

%%% make sure all paths end in a slash (forward vs. backwards depends on system)
if ispc==1
    sys_slash = '\';
else
    sys_slash = '/';
end
if work_path(end)~=sys_slash
    work_path = [work_path sys_slash];
end
if path_out_comp(end)~=sys_slash
    path_out_comp = [path_out_comp sys_slash];
end
if path_out_freq(end)~=sys_slash
    path_out_freq = [path_out_freq sys_slash];
end
for j=1:numel(path_raw_data)
    px = char(path_raw_data(j));
    if px(end)~=sys_slash
        px = [px sys_slash];
        path_raw_data(j) = cellstr(px);
    end
end