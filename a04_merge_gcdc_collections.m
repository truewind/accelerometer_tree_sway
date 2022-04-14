clear
clc
close all

%%% this script does one thing:
%
%%% (1) brings together sway frequency time series from multiple collection
%%% dates and merges into a single .mat file for each sensor/site

%% subset settings

%%% specify which sites you want to process (order defined in
%%% a00_data_file_settings.m
%%% can be single value or an array of numbers. count starts at 1.
sway_sites = [1 2];  % [1 2 3 4 5 6 7] 


%%% specify starting date and ending date for combined time series
sd_i = datenum(2014,10,1);  % yyyy,mm,dd
sd_f = datenum(2020,08,06);  % yyyy,mm,dd

dt = 3600;  % time step of time series (seconds) -- *** shoudl be same as fun_opt.win_int from a03_process_sway

%% load data and file settings

a00_data_file_settings;

%%% sub-select to only the sensors selected in sway_sites (above)
tree_names = tree_names(sway_sites);


%% ********** MAIN CODE **********

%%% gcdc fields
Gfields = {'Rstd', 'Dpts', 'fmax', 'fmax_pxx', 'f_pth', 'Amin', 'Amax', 'Amed', 'Aavg', 'Astd', 'Rmin', 'Rmax', 'Rmed', 'Ravg'};
ngf = numel(Gfields);


%%

dt = (dt/3600); % convert time step to hours

nTrees = numel(tree_names);

for j=1:nTrees

    %%% name of current sensor/tree
    xT = char(tree_names(j));
    
    %%% remove any previous version of the _all mat file
    allSwayFile = [path_out_freq 'GCDC_L02_Tree_Sway_' xT '_all.mat'];
    if exist(allSwayFile, 'file')==2
        delete(allSwayFile);
    end
    
    %%% find all .mat frequency analysis files for this sensor/tree.
    %%% specify as level 2 so we ignore level 1 (meta data) files, in case
    %%% stored in the same directory
    DD=dir([path_out_freq 'GCDC_L02*' char(xT) '*.mat']);

    %%% count the files
    nfiles = size(DD,1);
    
    %%% cycle through the files, load, and combine
    for k=1:nfiles
        ifile = DD(k).name;
        disp(ifile)
        
        %%% load GCDC data
        load([path_out_freq ifile])
        
        %%% initialize SWAY structure and variables during first loop
        if k==1
            SWAY.TIME = time_builder(sd_i, sd_f, dt);
            nt = size(SWAY.TIME,1);
            for qq=1:ngf
                eval(['SWAY.' char(Gfields(qq)) '=nan(nt,3);'])
            end
        end
        
        %%% grid data
        SWAY = data_grid(SWAY, GCDC, 1, 1, 0);
    end

    %%% save and clean up
    save(allSwayFile, 'SWAY')
    clear SWAY
end