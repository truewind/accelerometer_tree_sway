# accelerometer_tree_sway
Code used to process tree sway frequency from accelerometer data

## References
Raleigh, M.S., Gutmann, E.D., Van Stan II, J.T., Burns, S.P., Blanken, P.D., and E.E. Small (2022). Challenges and capabilities in estimating snow mass intercepted in conifer canopies with tree sway monitoring. *Water Resources Research*, doi: 10.1029/2021WR030972, https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021WR030972

## Data and Setup
* Acceleration data are collected with Gulf Coast Data Concepts (GCDC) accelerometers
* Code is written in Matlab R2020B. All files required to process the data are included in this repo.

## File structure for storing GCDC raw acceleration files
*Based on GCDC convention, the raw acceleration data are saved as .CSV files with name "DATA-XXX" where XXX= file number. These files are stored in separate folders, based on the collection date (i.e. when they were downloaded). Here the "collection date" is a subdirectory within the "path_raw_data" directory. These are named as YYYYMMDD where YYYY=year, MM=month, DD=date.*
  * path_raw_data
    * 20191218
      * DATA-001.CSV
      * DATA-002.CSV
      * DATA-003.CSV
      * ...
      * DATA-999.CSV
    * 20202012
      * DATA-001.CSV
      * DATA-002.CSV
      * DATA-003.CSV
      * ...
      * DATA-464.CSV 
    * 20200513
      * DATA-001.CSV
      * DATA-002.CSV
      * DATA-003.CSV
      * ...
      * DATA-751.CSV

## Steps
1. Configure paths and settings in **a00_data_file_settings.m**. If any paths do not exist, create them now.
2. Check the validity of the GCDC files and comile .CSV files into a single NetCDF for each collection:
    * open **a01_precheck_compile_GCDC.m**
    * change subset settings to choose your sensor/tree (order specified in **a00_data_file_settings.m**) and your date range
    * save and run **a01_precheck_compile_GCDC.m**
    * The code will throw an error if there are any strange dates (sometimes this is a quirk with GCDC files that is usually correctable). I usually correct the original raw .CSV file that has a bad date/time in the header by examining the files before/after and by keeping records on when I start/stop sensors.
    * If there are no errors, the script will produce a NetCDF of all acceleration data for each collection date in the time span specified.
3. Check the chronological order of the data in the NetCDF produced in Step 2, and sort if out of order:
    * open and run **a02_check_sorting.m**
    * messages will appear on screen as it inspects each NetCDF in the path where you saved the NetCDF files (specified in **a00_data_file_settings.m**)
4. Conduct a Lomb-Scargle frequency analysis of the acceleration data:
    * open **a03_process_sway.m"
    * change subset settings to choose your sensor/tree (order specified in **a00_data_file_settings.m**) and your date range
    * change frequency analysis settings and settings related to power specturm matrix visualization. Note! This includes the window size of the frequency analysis (e.g. a non-overlapping 5-minute window would require win_size=300 and win_int=300).
    * save and run **a01_precheck_compile_GCDC.m**
    * the script will run the frequency analysis and save the output to a .mat file in the directory you specified in the settings file (**a00_data_file_settings.m**)
5. Combine the time series of tree sway frequency from multiple collections into a single file:
    * open **a04_merge_gcdc_collections.m"
    * change subset settings to choose your sensor/tree (order specified in **a00_data_file_settings.m**) and your date range
    * set the output timestep of the time seres of tree sway (parameter dt in seconds)
    * save and run **a04_merge_gcdc_collections.m"
6. Filter and smooth the time series of tree sway frequency:
    * open **a05_filter_smooth_gcdc.m**
    * change subset settings to choose your sensor/tree (order specified in **a00_data_file_settings.m**) 
    * if needed, revisit the filtering and smooth parameters specified in teh settings file (**a00_data_file_settings.m**)
    * save and run **a05_filter_smooth_gcdc.m**
    * this will generate figures to show you the filtered and smoothed time series and the original data. you may consider iterating on this step and adjusting the filtering and smoothing parameters.
    * the script will output a .mat file with the smoothed and filtered time series that you can use for further analyses. The tree sway variables from various steps will be stored in a structure called "SWAY" with the following fields:
      * TIME = Lx7 matrix with time (1st col = year, 2nd col = month, 3rd col = day, 4th col = hour, 5th col = min, 6th col = julian date, 7th col = matlab serial date), where L = number of time steps   
      * fmax = Lx3 matrix with time series of tree sway frequency (Hz) for the three axes (1st col = X axis, 2nd col = Y axis, 3rd col = Z axis... see GCDC manual for orientations and your own notes on how the sensor was oriented in the field)
      * fmax_pxx = same as fmax, but only retaining points that met the required power level threshold
      * f_filt = same as fmax_pxx but after filtering
      * fsmooth = same as f_filt but after smoothing
      * favg = fsmooth averaged across axes (currently the Y and Z axis - modify as needed if different axes are of interest)
      * favg_gaps = same as favg but only during steps when valid data are present for all axes (modify which are selected, as needed)
