clear
clc
close all

%%% this script does two things:
%
%%% (1) takes merged datafile for each sensor/site and filters based on
%%%      power threshold and fmax and fmin 
%%% (2) smooths the data with a smoothing spline

%% subset settings

%%% specify which sites you want to process (order defined in
%%% a00_data_file_settings.m
%%% can be single value or an array of numbers. count starts at 1.
sway_sites = [1 2];  % [1 2 3 4 5 6 7] 


%%% other 
%% load data and file settings

a00_data_file_settings;

%%% sub-select to only the sensors selected in sway_sites (above)
tree_names = tree_names(sway_sites);
nTrees = numel(tree_names);


for j=1:nTrees
    
    %%% name of current sensor/tree
    xT = char(tree_names(j));
    
    %%% number of current tree
    iTree = sway_sites(j);
    fmax = fsway_max(iTree);
    fmin = fsway_min(iTree);
    
    %%% load file
    allSwayFile = [path_out_freq 'GCDC_L02_Tree_Sway_' xT '_all.mat'];
    load(allSwayFile);
    
    TIME=SWAY.TIME;
    
    
    %%% initialize variable for smooth and filt freq and gaps tracker
    SWAY.fsmooth = SWAY.fmax.*NaN;
    SWAY.f_filt  = SWAY.fmax.*NaN;
    SWAY.fgaps = SWAY.fmax.*NaN;
    
    for k=2:3
        %%% grab data for this axis
        f= SWAY.fmax(:,k);
        f_orig = f;
        fmax_pxx = SWAY.fmax_pxx(:,k);
        
        %%% remove points that do not meet the power level threshold plomb
        a = find(fmax_pxx>=pthresh);
        f(a) = NaN;

        %%% remove points that do not meet the max and min f
        f(f>fmax | f<fmin) = NaN;
        
        %%% remove any high outliers (3 day window, centered)
        TF = isoutlier(f,'movmean', [36 36]);
        idx_outlier = find(TF==1);
        f(idx_outlier) = NaN;
        
        %%% save the filtered data
        f_filt = f;
        
        %%% smooth the data and again apply max/min constrains
        [fitresult, gof] = createFit0(TIME(:,end),f, SmoothParam_sway);
        f_smooth = feval(fitresult, TIME(:,end));
        f_smooth(f_smooth>fmax | f_smooth<fmin) = NaN;
        
        
        %%% remove any spline results before first data or after last data
        firstdata = find(isnan(f)==0,1,'first');
        lastdata = find(isnan(f)==0,1,'last');
        f_smooth(1:firstdata-1) = NaN;
        f_smooth(lastdata+1:end) = NaN;
        

        figure
        hold on
        plot(f_orig, '.k')
        plot(idx_outlier, f_orig(idx_outlier), 'xk')
        plot(f, 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 3)
        plot(f_smooth)
        plot(zeros(size(f_orig))+fmax, '--k')
        plot(zeros(size(f_orig))+fmin, '--k')

        %%% store
        SWAY.f_filt(:,k) = f_filt;
        SWAY.fsmooth(:,k) = f_smooth;
        
        close all
        
        %%% keep track of missing values so we can break the spline at
        %%% places where the gaps are too big
        fgaps = isnan(f_filt);
        [fgaps_labeled, fgaps_num] = bwlabel(fgaps);
        fgaps_measure= regionprops(fgaps_labeled, fgaps, 'PixelValues');
        fgaps2 = fgaps_labeled;
        for fg=1:fgaps_num
            a=find(fgaps_labeled==fg);
            fgaps2(a) = numel(fgaps_measure(fg).PixelValues);
        end
        
        SWAY.fgaps(:,k) = fgaps2;
    end
    
    %%% averaged sway between axes
    SWAY.favg = nanmean(SWAY.fsmooth(:,2:3),2);
    
    fgaps = SWAY.fgaps;
    fgaps(SWAY.fgaps>gap_max) = NaN;
    fgaps(SWAY.fgaps<=gap_max) = 1;
    fgaps = mean(fgaps(:,2:3),2); % use mean, non nanmean so we only retain steps with both non-nan across axes
    SWAY.favg_gaps = SWAY.favg.*fgaps;
    
    %%% save file
    smoothSwayFile = [path_out_freq 'GCDC_L03_Tree_Sway_' xT '_smooth.mat'];
    save(smoothSwayFile, 'SWAY');
    
end