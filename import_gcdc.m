% imports Gulf Coast Data Concepts file into Matlab
%
% RELEASE NOTES
%   Written by Mark Raleigh (mraleig1@gmail.com), Feb 2018
%
% SYNTAX
%   GCDC = import_gcdc(GCDC_file)
%
% INPUTS
%   GCDC_file = str or cell string of file name (may include path) with .csv extension
%
% OUTPUTS
%   GCDC = structure with data, with fields:
%       SD = Lx1 array of serial dates of data
%       SC = Lx1 array of seconds elapsed
%       A0 = Lx1 array of acceleration on axis x (units = g's)
%       A1 = Lx1 array of acceleration on axis y (units = g's)
%       A2 = Lx1 array of acceleration on axis z (units = g's)
%       mTEMP = 1x1 value, metadata, temperature at start of file
%       mBATT = 1x1 value, metadata, battery at start of file
%       mSAMP = 1x1 value, metadata, sampling rate (Hz)
%       mSDi = 1x1 value, metadata, serial date of start of file
%

function GCDC = import_gcdc(GCDC_file)

GCDC_file = char(GCDC_file);

fid = fopen(GCDC_file);

zlock = 1;
while zlock>=1
    
    Q = fgetl(fid);
    if isempty(strfind(Q,'Start_time'))==0
        a = find(Q==',',1,'first');
        
        yr = str2double(Q(a+2:a+5));
        mo = str2double(Q(a+7:a+8));
        da = str2double(Q(a+10:a+11));
        HR = str2double(Q(a+14:a+15));
        MN = str2double(Q(a+17:a+18));
        SC = str2double(Q(a+20:a+25));
        SDi = datenum(yr,mo,da,HR,MN,SC);
        GCDC.mSDi = SDi;
        
        zlock=0;
        
        %%% get data
        M = importdata(GCDC_file);
        M = M.data;
        A0 = M(:,2)/2048;       % convert to g's (axis pointing down)
        A1 = M(:,3)/2048;       % convert to g's
        A2 = M(:,4)/2048;       % convert to g's

        SC = M(:,1);    % seconds
        SD = SDi+(M(:,1)/86400);
        
        if SDi>now
            % do this check because sometimes the GCDC
            % file saves a weird date (like July 2037)
            % - need to flag these files and rename
            % them so they are not included in the
            % processing
            disp('warning!! suspicious time step!')
        end
        
        GCDC.SD = SD;
        GCDC.SC = SC;
        GCDC.A0 = A0;
        GCDC.A1 = A1;
        GCDC.A2 = A2;
        
    else
        zlock = zlock+1;
    end
end
 
%%% continue reading metadata and extract any relevant variables. stop
%%% reading once we reach the data (updated by MSR, Feb 17, 2020)
zlock = 1;
while zlock>=1
    
    %     %%% then continue reading metadata
    %     Q = fgetl(fid);
    %     a = find(Q==',');
    %     GCDC.mTEMP = str2double(Q(a(1)+1:a(2)-1));
    %     GCDC.mBATT = str2double(Q(a(4)+1:a(5)-1));
    %
    %     Q = fgetl(fid);
    %     a = find(Q==',');
    %     GCDC.mSAMP = str2double(Q(a(1)+1:a(2)-1));
    %     GCDC.mSDi = SDi;
    Q = fgetl(fid);
    
    if numel(Q)>0
        if Q(1)==';'
            %%% then this is still a header line. see if there are any
            %%% varaibles of interest
            
            if isempty(strfind(Q,'Temperature'))==0
                a = find(Q==',');
                GCDC.mTEMP = str2double(Q(a(1)+1:a(2)-1));
            end
            
            if isempty(strfind(Q,'Vbat'))==0
                a = find(Q==',');
                GCDC.mBATT = str2double(Q(a(4)+1:a(5)-1));
            end
            
            if isempty(strfind(Q,'SampleRate'))==0
                a = find(Q==',');
                GCDC.mSAMP = str2double(Q(a(1)+1:a(2)-1));
            end
        else
            %%% then we are no longer on a metadata line. stop the while
            %%% loop
            zlock =0;
        end
    end
end

if exist('GCDC', 'var')~=1
    GCDC = [];
end

fclose(fid);