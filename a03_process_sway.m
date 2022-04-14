clear
clc
close all

%%% This script does one thing:
%
%%% (1) conducts a Lomb-Scargle frequency analyis of the acceleration data.
%%% settings are provided to specify a subset of sensors and collection
%%% dates, a subset of axes (GCDC records acceleration in 3 dimensions) to
%%% conduct the freq anaylsis, and other options related to Lomb-Scargle,
%%% windowing, and power spectrum visualization.


%% subset settings

%%% specify which sites you want to process (order defined in
%%% a00_data_file_settings.m
%%% can be single value or an array of numbers. count starts at 1.
sway_sites = 1; % [1 2];  % [1 2 3 4 5 6 7] 

% specify earliest and latest date of data collection (year, month, day),
% where data collection is when I downloaded all the data on the GCDC
% sensor and then cleared the memory on the sensor for the next period of
% data recording.
sd1 = datenum(2017,03,1); % earliest
sd2 = datenum(2017,03,12); % latest

%% frequency analysis settings

%%% indicate which of the three axes you want to do the frequency analysis
%%% enter 1 for axes you want to do the frequency analysis, 0 otherwise
%%% (1st var= Ax (usually down), 2nd= Ay, 3rd = Az)
f_axis = [0 1 1];

%%% enter the sampling frequency (Hz) at which the sensor sampled
SAMPFREQ = 12;

%%% Lomb-Scargle periodogram parameters (see plomb.m documentation)
fun_opt.ofac = 2;   % oversampling factor
fun_opt.fhi = 3;    % maximum frequency (Hz) to consider
fun_opt.Pd = 0.99;  % power level threshold (probability of a true signal peak)

%%% windowing settings
%%% The frequency analysis will be conducted in independent windows. All
%%% available data within the window will be used in the Lomb-Scargle
%%% frequency analysis.
fun_opt.win_size=300;    % window width (seconds)
fun_opt.win_int=300;     % window spacing (seconds)


%%% note if win_int<win_size, then you will get overlapping windows (i.e.,
%%% data that is used in multiple windows) and possibly a smoother
%%% response, but will incur more computation time than win_size=win_int.


%% settings related to power spectrum matrix (2d with dimensions of time and frequency bin)

%%% create vector for power spectrum visualization
pinc = 0.01;    % increment in frequency value (Hz)
pmin = 0.4;     % minimum frequency value (Hz)
pmax = 6.0;     % maximum frequency value (Hz)
fun_opt.fhi = 6;    % maximum frequency (Hz) to consider

%%% build the vector
pvec = pmin-(pinc/2):pinc:pmax+(pinc/2);
pvec = [0-(pinc/2), pvec, fun_opt.fhi + (pinc/2)];
np = numel(pvec)-1;

%% load data and file settings

a00_data_file_settings;

%%% sub-select to only the sensors selected in sway_sites (above)
path_raw_data = path_raw_data(sway_sites);
tree_names = tree_names(sway_sites);




%% ********** MAIN CODE **********
%%% (in general, should not need to make any changes below)


%%% constants
gconv = 2048; % acceleration conversion factor from GCDC manual
secondsday = 86400; % seconds in a day

%%% find .nc files in data path
dir_rec = dir(path_out_comp);
dir_rec = dir_rec(3:end);
nc_files = {''};
nnc = 0;

sens_subset = [];
file_subset = {};
for j=1:size(dir_rec,1)
    a = strfind(dir_rec(j).name, '.nc');
    if isempty(a)==0
        nnc=nnc+1;
        nc_files(nnc,1) = cellstr(dir_rec(j).name);
        
        nc_currFile = char(nc_files(nnc));

        %%% check if this file is within the subset (first name, then date)
        iFlag = 0;
        for k=1:numel(tree_names)
            if ~contains(nc_currFile, char(tree_names(k)))==0
                iFlag = iFlag+1;
            end
        end
        if iFlag==1
            %%% get collection date of this .nc file
            cYYYY = str2double(nc_currFile(end-10:end-7));
            cMM = str2double(nc_currFile(end-6:end-5));
            cDD = str2double(nc_currFile(end-4:end-3));
            cDATE = datenum(cYYYY,cMM,cDD);
            if cDATE>=sd1 && cDATE<=sd2+1
                iFlag = iFlag+1;
            end
        end
        
        if iFlag==2
            %%% then this sensor is listed in the subset settings
            sens_subset = [sens_subset, j];
            file_subset(numel(sens_subset),1) = cellstr(nc_currFile);
        end

    end
end

disp(' ')
disp(' ')
disp('Subset Settings:')
disp('*** Sensors/Trees:')
disp(tree_names)
disp(' ')
disp('*** Date of collection:')
disp(['Between ' datestr(sd1,'yyyy-mm-dd') ' and ' datestr(sd2, 'yyyy-mm-dd')])
disp(' ')

if isempty(sens_subset)==1
    error('None of the sensors meet the subset settings. Exiting.') 
else
    disp(['.... found ' num2str(numel(file_subset)) ' collections meeting the subset settings'])
    disp(' ')
    disp('The following collections will be processed:')
    disp(file_subset)
    disp(' ')
    disp(' ')
    disp('Press any key to continue')
    pause
    
    nc_files = file_subset;
    nnc=numel(nc_files);
end

tic;

%%% cycle through nc files and conduct frequency analysis
for j=1:nnc
    
    %%% current nc file
    ncx = char(nc_files(j));
    nc_filepath = [path_out_comp ncx];
    
    %%% get serial dates and Accel data from .nc file
    ncSD = ncread(nc_filepath, 'serial_date');
    
    %%% find range in time for the data in this collection
    sdmin = nanmin(ncSD);
    [mY1,mM1,mD1,mH1,~,~] = datevec(sdmin);
    sdmin = datenum(mY1,mM1,mD1,mH1,0,0);
    
    nci=ncinfo(nc_filepath, 'serial_date');
    sdmax = nanmax(ncSD);
    [mY2,mM2,mD2,mH2,~,~] = datevec(sdmax);
    sdmax = datenum(mY2,mM2,mD2,mH2,0,0);
    
    %%% start new structure and initialize outputs
    GCDC = [];
    GCDC.TIME = time_builder(sdmin, sdmax, (fun_opt.win_int/3600) );
    nz = size(GCDC.TIME,1);
    GCDC.Dpts = nan(nz,3);
    GCDC.fmax = nan(nz,3);
    GCDC.fmax_pxx = nan(nz,3);
    GCDC.f_pth = nan(nz,3);
    
    %%% stats on detrended, demeaned accel data over each window
    GCDC.Amin = nan(nz,3);
    GCDC.Amax = nan(nz,3);
    GCDC.Amed = nan(nz,3);
    GCDC.Aavg = nan(nz,3);
    GCDC.Astd = nan(nz,3);
    
    %%% stats on RAW accel data (no detrend or demean) over each window
    GCDC.Rmin = nan(nz,3);
    GCDC.Rmax = nan(nz,3);
    GCDC.Rmed = nan(nz,3);
    GCDC.Ravg = nan(nz,3);
    GCDC.Rstd = nan(nz,3);
    
    %%% variables related to power spectrums
    GCDC.pxx_fcen = pvec(1:end-1)+(pinc/2);
    GCDC.pxx = nan(nz,np,2);
    
    %%% cycle through time steps
    for n=1:nz
        disp([ncx ': processing step ' num2str(n) '/' num2str(nz)])
        
        %%% serial dates of window limits
        w1 = GCDC.TIME(n,7);
        w2 = w1 - (fun_opt.win_size/86400);
        
        %%% find all files that have data for this timestep
        fff = find( ncSD >=w2 & ncSD <w1);
        
        %%% only continue if we have data in this range
        if isempty(fff)==0
            
            fstart = nanmin(fff);
            fcount = numel(fff);

            ncA0 = ncread(nc_filepath, 'Ax',fstart,fcount,1);
            ncA1 = ncread(nc_filepath, 'Ay',fstart,fcount,1);
            ncA2 = ncread(nc_filepath, 'Az',fstart,fcount,1);
            
            %%% initialize
            CDATA.SC = (ncSD(fff) - nanmin(ncSD(fff))) .* secondsday;  % seconds since start of period
            CDATA.A0 = double(ncA0);  % accel 0
            CDATA.A1 = double(ncA1);  % accel 1
            CDATA.A2 = double(ncA2);  % accel 2
            lombdata = [];
            
            %%% convert to G's
            CDATA.A0 = CDATA.A0 ./ gconv;
            CDATA.A1 = CDATA.A1 ./ gconv;
            CDATA.A2 = CDATA.A2 ./ gconv;
            
            %%% sort
            [~,iX] = sort(CDATA.SC);
            CDATA.SC = CDATA.SC(iX,:);
            CDATA.A0 = CDATA.A0(iX,:);
            CDATA.A1 = CDATA.A1(iX,:);
            CDATA.A2 = CDATA.A2(iX,:);
            
            %%% cycle through axes and analyze
            for fa=1:numel(f_axis)
                
                disp(['... processing axis ' num2str(fa)])
                
                if fa==1
                    lombdata = [CDATA.SC, CDATA.A0];
                elseif fa==2
                    lombdata = [CDATA.SC, CDATA.A1];
                elseif fa==3
                    lombdata = [CDATA.SC, CDATA.A2];
                end
                
                
                %%% get rid of any NaN rows (in time or data)
                a = find(isnan(nansum(lombdata,2))==0);
                lombdata = lombdata(a,:);
                
                %%% make sure it is sorted in chronological order
                lombdata=sortrows(lombdata,1);
                
                %%% eliminate duplicates
                [~,ls_uni,~]=unique(lombdata(:,1),'rows');
                lombdata = lombdata(ls_uni,:);
                
                %%% narrow down to only this timestep
                a = find(lombdata(:,1)>=0 & lombdata(:,1)<=fun_opt.win_size);
                lombdata=lombdata(a,:);
                
                %%% subtract mean and detrend
                pl_x = lombdata(:,2);
                raw_x = pl_x;
                pl_x = pl_x-nanmean(pl_x); % subtract the mean before LS
                pl_t = lombdata(:,1);
                pl_x = detrend(pl_x);
                
                %%% frequency analysis
                if f_axis(fa)==1
                    if isempty(lombdata)==1
                        error('check this')
                    else
                        %%% lomb-scargle
                        [ls_pxx,ls_f,pth]=plomb(pl_x,pl_t,fun_opt.fhi,fun_opt.ofac, 'Pd', fun_opt.Pd);
%                         error('fvec MUST be called fvec apparently if it is used as an input')
%                         hh = hann(numel(pl_x)); [ls_pxx,ls_f,pth]=plomb(pl_x.*hh,pl_t,6,fun_opt.ofac, 'Pd', fun_opt.Pd);
%                         fvec=logspace(-2,1,10); [ls_pxx,ls_f]=plomb(pl_x,pl_t,fvec);
%                         fvec=linspace(0.01,6,10); [ls_pxx,ls_f]=plomb(pl_x,pl_t,fvec);
                        
                        
%                         [ls_pxx,ls_f,pth]=plomb(pl_x,pl_t,6,fun_opt.ofac, 'Pd', fun_opt.Pd);
%                         [ls_pxx,ls_f]=plomb(pl_x,pl_t,0.3:0.01:6);
%                         [ls_pxx,ls_f]=plomb(pl_x,pl_t,0.01:0.01:6);
%                         ls_pxx = pow2db(ls_pxx);
                        
                        %%% make sure the output is sorted (but already done, but just in case)
                        [~,ls_I] = sort(ls_f);
                        ls_pxx = ls_pxx(ls_I);
                        ls_f = ls_f(ls_I);
                        
                        %%% apply filter (moving window)
                        ls_pxx=movingmean(ls_pxx,13,1,1);   % doesn't seem to make a difference
%                         ls_pxx=smoothdata(ls_pxx,'gaussian',[15 15], 'omitnan');   
%                         ls_pxx(ls_f<0.33) = NaN;
                        
                        %%% find max freq
                        ls_pxx2 = ls_pxx;
                        ls_pxx2(ls_f<0.33) = NaN;
                        mx_pxx = nanmax(ls_pxx2);
%                         mx_pxx = nanmax(ls_pxx);
                        fq_pxx = find(ls_pxx==mx_pxx);
                        GCDC.fmax(n,fa) = nanmean(ls_f(fq_pxx));
                        GCDC.fmax_pxx(n,fa) = mx_pxx;
                        GCDC.f_pth(n,fa) = pth;
                        
                        GCDC.Dpts(n,fa) = numel(pl_x);
                        
                        %%%
%                         figure(1);
% %                         hold on;
% %                         hold on;
%                         plot(ls_f, ls_pxx, '-k')
% % % % %                         plot(ls_f(fq_pxx), pDB(fq_pxx), 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
%                         set(gca, 'xscale', 'log')
% %                         set(gca, 'yscale', 'log')
% % % % %                         hold off;
%                         xlim([0.01 6])
% % % %                         ylim([-80 -35])
%                         xlabel('frequency [Hz]')
%                         ylabel('Py(f) [dB]')
%                         title(datestr(w1, 'yyyy-mm-dd HH:MM:SS'))
%                         drawnow;

                        
                        %%% store power spectrum for visualization
                        [~,~,ibin] = histcounts(ls_f,pvec);
                        
                        % enforce start/end limits on truncated hist
                        if ibin(1) ==0
                            ibin2 = ibin;
                            ibin2(ibin2<=0) = NaN;
                            ibin(1) = nanmin(ibin2);
                        end
                        
                        if ibin(end)==0
                            ibin2 = ibin;
                            ibin2(ibin2<=0) = NaN;
                            ibin(end) = nanmax(ibin2);
                        end
                        
                        valh = accumarray(ibin, ls_pxx, [numel(pvec)-1 1], @nanmean);
                        GCDC.pxx(n,:,fa) = valh';
                    end
                    
                    
                end
                
                %%% summarize accel. data
                GCDC.Amin(n,fa) = nanmin(pl_x);
                GCDC.Amax(n,fa) = nanmax(pl_x);
                GCDC.Amed(n,fa) = nanmedian(pl_x);
                GCDC.Aavg(n,fa) = nanmean(pl_x);
                GCDC.Astd(n,fa) = nanstd(pl_x);
                
                GCDC.Rmin(n,fa) = nanmin(raw_x);
                GCDC.Rmax(n,fa) = nanmax(raw_x);
                GCDC.Rmed(n,fa) = nanmedian(raw_x);
                GCDC.Ravg(n,fa) = nanmean(raw_x);
                GCDC.Rstd(n,fa) = nanstd(raw_x);
            end
        end
    end
    
    %%% saving output
    GCDC.meta_proc = fun_opt;
    
    matname = ncx; %'GCDC_L01_Raw_Data_GM08_collection_20170226.nc';
    a1 = find(matname == '_');
    a1 = a1(4)+1;
    a2 = find(matname == '.');
    a2 = a2(end)-1;
    matname = ['GCDC_L02_Tree_Sway_' matname(a1:a2) '.mat'];
%     matname = ['GCDC_L02_Tree_Sway_' matname(a1:a2) '_5min.mat'];
    save([path_out_freq matname], 'GCDC')
end


toc;




