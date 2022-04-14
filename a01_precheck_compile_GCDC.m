clear
clc
close all

%%% This script does two things:
%
%%% (1) checks for any GCDC data files that have a bad date. In my
%%% experience this typically happens randomly where a single file will
%%% suddenly be in the distant future (like 2037). The file before and
%%% after have the correct date, so it is usually possible to correct the wrong
%%% date. I do that correction manually, and this code tells me when it
%%% finds files that need to be corrected.
%
%%% (2) For every collection date, it compiles all raw GCDC data (which can
%%% be hundreds of .csv files) into a single netCDF file. No frequency
%%% analysis is done yet (that is done in a03_process_sway.m).

%% subset settings

%%% specify which sites you want to check (order defined in
%%% a00_data_file_settings.m
%%% can be single value or an array of numbers. count starts at 1.
sway_sites = 2;  % [1 2 3 4 5 6 7] 

% specify earliest and latest date of data collection (year, month, day),
% where data collection is when I downloaded all the data on the GCDC
% sensor and then cleared the memory on the sensor for the next period of
% data recording.
sd1 = datenum(2018,05,13); % earliest
sd2 = datenum(2018,05,17); % latest


%% load data and file settings

a00_data_file_settings;

%% ********** MAIN CODE **********
%%% (in general, should not need to make any changes below)

%%% sub-select to only the sensors selected in sway_sites (above)
path_raw_data = path_raw_data(sway_sites);
tree_names = tree_names(sway_sites);

%%% change initial path (necessary for running in the background?)
cd(work_path);

%% precheck for any weird data (when GCDC have a glitch, it is a single file that has the wrong time stamp in the future)

nsensors = numel(tree_names);

disp(' **** Checking for suspicious timestamps **** ')
for s=1:nsensors
    
    %%% tree name
    Tname = char(tree_names(s));
    
    %%% path to data of current sensor
    dpath = char(path_raw_data(s));
    
    %%% get listing of sub folders
    dir_rec = dir(dpath);
    dir_rec = dir_rec(3:end,:);
    
    
    %%% cycle through sub folders
    for m=1:size(dir_rec,1)
        %%% folder name
        FNAME = dir_rec(m,1).name;
        
        disp(['checking ' Tname ', in folder ' FNAME '...'])
        
        dir_rec2 = dir([dpath FNAME '/']);
        
        yrX = str2double(FNAME(1:4));
        
        SDmax = nanmax([dir_rec2.datenum].');
        SDmin = nanmin([dir_rec2.datenum].');
        
        if SDmax>now || SDmin < yrX-1
            %%% then check to make sure we have corrected the time header
            if SDmax>now
                a = find([dir_rec2.datenum].'==SDmax);
            else
                a = find([dir_rec2.datenum].'==SDmin);
            end
            
            a = a+1; % usually the file preceeding the bad file has the incorrect date modified/created
            
            
            GCDCtest = import_gcdc([dpath FNAME '/' char(dir_rec2(a).name)]);
            
            if GCDCtest.mSDi>now || GCDCtest.mSDi < yrX-1
                disp([dpath FNAME '/' char(dir_rec2(a).name)])
                error('... detected bad timestep at/near the above file!')
            else
                % assumed it was corrected manually already
                disp([dpath FNAME '/ contains files with corrected time headers - these passed'])
            end
        else
            disp('... passed')
        end
        
    end
end

disp('   ')
disp('   ')
disp('   ')
disp('   ')
disp('   ')

%% compile

disp(' **** Compiling Raw Data **** ')

for s=1:nsensors
    
    %%% tree name
    Tname = char(tree_names(s));
    
    %%%
    dpath = char(path_raw_data(s));

    %%% get listing of sub folders
    dir_rec = dir(dpath);
    
    
    
    %%% cycle through sub folders
    for m=1:size(dir_rec,1)
        
        %%% folder name
        FNAME = dir_rec(m,1).name;
        FFOLD = [dpath char(FNAME) '/'];
        
        
        
        %%% manual corrections to apply
        SD_Shift = 0;   % assume no shift is needed in time (some cases where the time step was off by an hour)
        
        if strcmp('Niwot_Tree02', Tname)==1 && strcmp(char(FNAME), '20150101')==1
            SD_Shift = datenum(2014,11,6)-datenum(2014,10,6);   % accidentally started on wrong month (off my 1). move forward
        elseif strcmp(char(FNAME), '20161125')==1
            if strcmp('Niwot_Tree01', Tname)==1
                SD_Shift = ((24/60)+15)/24; % 15 hrs , 24 min behind approximately
            elseif strcmp('Niwot_Tree02', Tname)==1
                SD_Shift = ((23/60)+15)/24; % 15 hrs , 23 min behind approximately
            end
        else
            SD_Shift = 0;
        end
        
        
        
        %%% check if folder name is expected length (8 characters)
        if numel(char(FNAME))==8
            %%% get subdir listing, filter any 0 size files out
            subdir_rec=dir(FFOLD);
            subdir_rec = subdir_rec(3:end);
            sub_files_size = [subdir_rec.bytes].';
            a = find(sub_files_size>0);
            subdir_rec = subdir_rec(a);
            nfiles = size(subdir_rec,1);
            
            SD_collection = datenum(str2double(FNAME(1:4)), str2double(FNAME(5:6)), str2double(FNAME(7:8)));
            
            %%% check if the collection date is in range and if the folder has files
            if SD_collection>=sd1 && SD_collection<=sd2 && nfiles>=1
                
                %%% initialize metadata structure
                GCDC_meta.SD = [];
                GCDC_meta.mTEMP = [];
                GCDC_meta.mBATT = [];
                GCDC_meta.mSAMP = [];
                
      
                %%% counter
                Z=0;
                
                coll_name = ['GCDC_L01_Raw_Data_' char(tree_names(s)) '_collection_' char(FNAME)];
                

                %%% NETCDF initialize
                if Z==0
                    nc_file = [path_out_comp coll_name '.nc'];
                    if exist(nc_file, 'file')==2
                        delete(nc_file);
                    end
                    
                    ncid = netcdf.create(nc_file, 'NC_WRITE');
                    dimid = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));
                    
                    varid1 = netcdf.defVar(ncid, 'serial_date', 'NC_DOUBLE', dimid); 

                    varid2 = netcdf.defVar(ncid, 'Ax', 'NC_SHORT', dimid);
                    varid3 = netcdf.defVar(ncid, 'Ay', 'NC_SHORT', dimid);
                    varid4 = netcdf.defVar(ncid, 'Az', 'NC_SHORT', dimid);
                    
                    netcdf.endDef(ncid);
                    
                end
                

                
                for n=1:nfiles
                    disp(['Now importing ' num2str(n) '/' num2str(nfiles) ' in folder ' char(FNAME) ' for tree ' Tname])
                    
                    %%% current file number
                    FNi = subdir_rec(n).name;
                    GFILE = [FFOLD FNi];
                    
                    %%% only open if 12 char long file name (DATA-XXX.CSV)
                    if numel(FNi)==12
                        %%% open file
                        M = import_gcdc(GFILE);
                        Mn = size(M.SD,1);
                        
                        %%% apply any shift (correction)
                        if SD_Shift~=0
                            disp(['.... applying manual time correction of ' num2str(SD_Shift*24) ' hours'])
                            M.SD = M.SD + SD_Shift;
                            M.mSDi = M.mSDi + SD_Shift;
                        end
                        
                        %%% store metadata
                        GCDC_meta.SD = [GCDC_meta.SD; M.mSDi];
                        GCDC_meta.mTEMP = [GCDC_meta.mTEMP; M.mTEMP];
                        GCDC_meta.mBATT = [GCDC_meta.mBATT; M.mBATT];
                        GCDC_meta.mSAMP = [GCDC_meta.mSAMP; M.mSAMP];
                        
                        
                        %%% time
                        time = M.SD;
                        
                        %%% convert back to voltage integers for lower storage
                        Ax = M.A0.*2048;
                        Ay = M.A1.*2048;
                        Az = M.A2.*2048;
                        
                        Ax = int16(Ax);
                        Ay = int16(Ay);
                        Az = int16(Az);
                        
                        Ax(isnan(Ax)) = -99999;
                        Ay(isnan(Ay)) = -99999;
                        Az(isnan(Az)) = -99999;
                        
                        if Z==0
                            netcdf.putVar(ncid, varid1, 0, Mn, M.SD);
                            netcdf.putVar(ncid, varid2, 0, Mn, Ax);
                            netcdf.putVar(ncid, varid3, 0, Mn, Ay);
                            netcdf.putVar(ncid, varid4, 0, Mn, Az);
                            
                        else
                            %%% now create a netcdf with the raw data
                            ncid = netcdf.open(nc_file, 'NC_WRITE');
                            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,0);
                            
                            varid = netcdf.inqVarID(ncid, 'serial_date');
                            [dimname, dimlen] = netcdf.inqDim(ncid,0);
                            netcdf.putVar(ncid, varid, dimlen, Mn, time);
                            
                            varid = netcdf.inqVarID(ncid, 'Ax');
                            netcdf.putVar(ncid, varid, dimlen, Mn, Ax);
                            
                            varid = netcdf.inqVarID(ncid, 'Ay');
                            netcdf.putVar(ncid, varid, dimlen, Mn, Ay);
                            
                            varid = netcdf.inqVarID(ncid, 'Az');
                            netcdf.putVar(ncid, varid, dimlen, Mn, Az);
                        end
                        
                        
                        netcdf.close(ncid);
                        Z = Z+1;
                        
                        
                        
                    end
                    
                end
                
                %%% save metadata
                save([path_out_comp coll_name '_meta.mat'], 'GCDC_meta');
            end
        end
        
    end
    
end



