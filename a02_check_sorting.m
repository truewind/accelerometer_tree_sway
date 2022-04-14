clear
clc
close all

%%% This script does one thing:
%
%%% (1) checks the compiled netCDF files that were produced with
%%% a01_precheck_compile_GCDC.m to ensure that each is sorted
%%% chronologically. If there are places where the data are out of
%%% chronological order, it will sort the data and save it (overwriting the
%%% netCDF to a new one).

%% load data and file settings

a00_data_file_settings;



%% ********** MAIN CODE **********
%%% (in general, should not need to make any changes below)


%%

D = dir(path_out_comp);
D = D(3:end);

n = size(D,1);


for s=1:n
    
    %%%
    FX = char(D(s).name);
    
    
    if contains(FX, '.nc')==1
        
        %%% load serial date of nc file
        nc_file = [path_out_comp FX];
        serial_date = ncread(nc_file, 'serial_date');
        sd_min = nanmin(serial_date);
        sd_max = nanmax(serial_date);
        
        TF = issorted(serial_date);
        
        disp(nc_file);
        if TF==0
            disp('... not sorted!')
            
            %%% find sorting indices
            [~,I]=sort(serial_date);
            
            %%% load Accel data
            Ax = ncread(nc_file, 'Ax');
            Ay = ncread(nc_file, 'Ay');
            Az = ncread(nc_file, 'Az');
            
            %%% sort
            serial_date = serial_date(I);
            Ax = Ax(I);
            Ay = Ay(I);
            Az = Az(I);
            clear I
            
            %%% check if sorted again
            TF2 = issorted(serial_date);
            
            if TF2==0
                error('.... failed to sort the .nc file')
            else
                %%% write new nc_file
                delete(nc_file);

                ncid = netcdf.create(nc_file, 'NC_WRITE');
                dimid = netcdf.defDim(ncid, 'time', netcdf.getConstant('NC_UNLIMITED'));
                
                varid1 = netcdf.defVar(ncid, 'serial_date', 'NC_DOUBLE', dimid);
                varid2 = netcdf.defVar(ncid, 'Ax', 'NC_SHORT', dimid);
                varid3 = netcdf.defVar(ncid, 'Ay', 'NC_SHORT', dimid);
                varid4 = netcdf.defVar(ncid, 'Az', 'NC_SHORT', dimid);               
                netcdf.endDef(ncid);
                
                Mn = size(serial_date,1);
                netcdf.putVar(ncid, varid1, 0, Mn, serial_date);
                netcdf.putVar(ncid, varid2, 0, Mn, Ax);
                netcdf.putVar(ncid, varid3, 0, Mn, Ay);
                netcdf.putVar(ncid, varid4, 0, Mn, Az);
                
                netcdf.close(ncid);
                
                clear serial_date Ax Ay Az
                disp('---------- sorted now! -------------')
                disp(' ')
            end

        else
            disp('... already sorted')
        end
        
        disp([      'first date=' datestr(sd_min, 'yyyy-mm-dd HH:MM:SS')])
        disp([      'last date=' datestr(sd_max, 'yyyy-mm-dd HH:MM:SS')])
    end
    
    
end






