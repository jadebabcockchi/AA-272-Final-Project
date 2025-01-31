function [gps_ephem, ionoGAL] = read_clean_GALbroadcast(navfilename,clean)

%==========================================================================
%==========================================================================
% [gps_ephem,ionoparams] = read_GPSbroadcast(navfilename,clean)
%
% Reads an IGS GPS Broadcast Ephemeris for all GPS satellites and 
%  constructs a matrix of all ephemeris values. Note that each row of the 
%  gps_ephem is a new entry in the broadcast file which will most likely 
%  contain multiple entries for the same PRN.
%
%
% Author: Ben K. Bradley
% Date: 07/19/2009
%
%
% INPUT:             Description                                      Units
%
%  navfilename   - name of IGS broadcast ephemeris file to read in   string
%  clean - optional flag to clean the incorrect entries from the ephem 
%          default is false.  If you want to clean, call it with true.
%
% OUTPUT:       
%    
%  gps_ephem     - matrix of gps satellite orbit parameters          (nx25)
%  
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: delta_n, mean motion difference from computed value, rad/s
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: i_rate, rate of change of inclination angle, rad/s
%                 col11: Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
%                 col12: Cus, amplitude of the sine harmonic correction term to the argument of latitude
%                 col13: Crc, amplitude of the cosine harmonic correction term to the orbit radius
%                 col14: Crs, amplitude of the sine harmonic correction term to the orbit radius
%                 col15: Cic, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col16: Cis, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col17: Toe, reference time ephemeris (seconds into GPS week)
%                 col18: IODE, issue of data (ephemeris) 
%                 col19: GPS_week, GPS Week Number (to go with Toe) wrt 06jan80 (no rollover)
%                 col20: Toc, time of clock
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: Af2, satellite clock drift rate (sec/sec/sec)
%                 col24: TGD, satellite group delay to be used for single f
%                 col25: health, satellite health (0=good and usable)
%
%
%  ionoparams    - parameters for the Klobuchar  [A0 A1 A2 A3 B0 B1 B2 B3]
%                    ionospheric model                        
%          
%
% Coupling:
%
%  none
%
% References:
% 
%  [1] Interface Control Document: IS-GPS-200M - page 105
%        < https://www.gps.gov/technical/icwg/IS-GPS-200M.pdf >
%
%  [2] RINEX GPS Format, Version 2, (Table A4) - page A23
%        < https://files.igs.org/pub/data/format/rinex304.pdf >
%
%==========================================================================
%==========================================================================


%% Header

% Open desired IGS Ephemeris File
%==========================================================================

if (exist(navfilename,'file') == 2)
    fid = fopen(navfilename,'r');
else
    error(sprintf('Unable to find broadcast file: %s',navfilename), 'ERROR!');
end

if nargin < 2
    clean = false;
end

% Step through the header of the file and pull out iono parameters
%==========================================================================
headerend   = [];
ALPHA = [];
BETA  = [];

while (isempty(headerend) == 1)
   tline = fgetl(fid); 
   headerend = findstr(tline, 'END OF HEADER');
   
   % Check for IONOSPHERIC CORR lines
   ionoGPSA = findstr(tline, 'GPSA');

   % Check for IONOSPHERIC CORR lines
   ionoGPSB = findstr(tline, 'GPSB');
   
   if (isempty(ionoGPSA) == 0)
       % The line contains ionospheric correction data
       cellA = strsplit(tline);
       A0 = cell2mat(cellA(2));
       A1 = cell2mat(cellA(3));
       A2 = cell2mat(cellA(4));
       A3 = cell2mat(cellA(5));
       
       % Store alpha parameters for this satellite
       ALPHA = [str2num(A0) str2num(A1) str2num(A2) str2num(A3)];
   end
   
   % Similarly check for IONOSPHERIC CORR for the BETA parameters
   if (isempty(ionoGPSB) == 0)
       % The line contains ionospheric correction data
       cellB = strsplit(tline);
       B0 = cell2mat(cellB(2));
       B1 = cell2mat(cellB(3));
       B2 = cell2mat(cellB(4));
       B3 = cell2mat(cellB(5));
       
       % Store alpha parameters for this satellite
       BETA = [str2num(B0) str2num(B1) str2num(B2) str2num(B3)];
       
       % Store the satellite name and BETA parameters
       ionoGAL = [ALPHA BETA];
   end
end


%% Body

j = 1;

% Enter main loop to read the rest of the ephemeris file
%==========================================================================
%==========================================================================
while 1
    
    % Load next line in ephemeris file
    tline = fgetl(fid);
    
    
    % If the next line is not a character then the end of the file has been
    %   reached and the while loop is exited
    if ~ischar(tline)
        break; 
    end

    if tline(1) == 'E'
   

        %-----------------------------------------------------------------
        % Read in variables of the FIRST line of this satellite's ephemeris
        %-----------------------------------------------------------------
        prn = str2num(tline(2:3));   % PRN number of satellite
       
        Af0 = str2num(tline(25:42)); % clock bias  (s)
        Af1 = str2num(tline(44:61)); % clock drift (s/s)
        Af2 = str2num(tline(63:80)); % clock drift rate (s/s^2)
        

        
        %-----------------------------------------------------------------
        % SECOND LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in second line of satellite ephemeris
        
        IODE = str2num(tline(6:23));    % Issue of Data (Ephemeris)  
        
        Crs  = str2num(tline(25:42));   % Amplitude of the Sine Harmonic Correction 
                                        %  Term to the orbit radius, meters
      
        delta_n= str2num(tline(44:61)); % Mean Motion Difference from Computed Value, rad/s
        
        M0   = str2num(tline(63:80));   % Mean Anomaly at Reference Time, rad
        
        %-----------------------------------------------------------------
        % THIRD LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in third line of satellite ephemeris
        
        Cuc = str2num(tline(6:23));     % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Argument of
                                        %  Latitude, rad
       
        ecc = str2num(tline(25:42));    % Eccentricity
        
        Cus = str2num(tline(44:61));    % Amplitude of the Sine Harmonic Correction
                                        %  Term to the Argument of
                                        %  Latitude, rad
       
        sqrt_a = str2num(tline(63:80)); % Square root of the semi-major Axis, sqrt(meters)
        
        %-----------------------------------------------------------------
        % FOURTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in fourth line of satellite ephemeris
        
        Toe = str2num(tline(6:23));     % Reference Time Ephemeris (sec into GPS week)
        
        Cic = str2num(tline(25:42));    % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Angle of
                                        %  Inclination, rad
        
        Loa = str2num(tline(44:61));    % Longitude of Ascending Node of Orbit Plane
                                        %  at Weekly Epoch, rad
      
        Cis = str2num(tline(63:80));    % Amplitude of the Sine Harmonic Correction
                                        %  Term to the Angle of
                                        %  Inclination, rad
         
        %-----------------------------------------------------------------                    
        % FIFTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in fifth line of satellite ephemeris
        
        incl = str2num(tline(6:23));    % Inclination Angle at Reference Time, rad
        
        Crc  = str2num(tline(25:42));   % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Orbit Radius, meters
     
        perigee = str2num(tline(44:61));% Argument of Perigee, rad
        
        ra_rate = str2num(tline(63:80));% Rate of Change of Right Ascension, rad/s
        
        %-----------------------------------------------------------------
        % SIXTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in sixth line of satellite ephemeris
        
        i_rate = str2num(tline(6:23));   % Rate of change of inclination angle, rad/s
        
        %str = tline(23:41);             % codes on L2 channel (unecessary)
       
        GPS_week = str2num(tline(44:61));% GPS Week Number (to go with Toe)
        if GPS_week < 1024
            GPS_week = GPS_week+2048; % Penny added to avoid modulo 1024
        end
        %str   = tline(61:79);           % L2 flag
        
        %-----------------------------------------------------------------
        % SEVENTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Includes: SV accuracy, SV health, TGD, IODC
        
        health = str2num(tline(24:42)); % Satellite health (0.00 = usable)
        TGD    = str2num(tline(44:61)); % TGD Timing Group Delay, seconds

        
        %-----------------------------------------------------------------
        % EIGHTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in eighth line of satellite ephemeris
        
        Toc = Toe; % Time of clock
        
        
        
        
        gps_ephem(j,:) = [prn M0 delta_n ecc sqrt_a Loa incl perigee ra_rate i_rate Cuc Cus Crc Crs Cic Cis Toe IODE GPS_week Toc Af0 Af1 Af2 TGD health];
        
        j = j + 1;  
    end
        
end

if clean
    prn_all = unique(gps_ephem(:,1));  % Need to check them by PRN
    for iprn = 1:length(prn_all)
        prn = prn_all(iprn);
    % Pull out ephemerides for PRN in question
    kk  = find(gps_ephem(:,1) == prn);  % kk is vector containing row numbers of ephem_all that are for sat.no. 'index' 
    sat_ephem = gps_ephem(kk,:);        % sat_ephem is matrix of all ephem data for each entry of sat.no. 'index'
    % No matching PRN found, returning data will be NaNs
    if ~isempty(kk) 
    % Remove bad ephemerides. Occasionally, new ephemerides are uploaded at odd
    %  times which occur prior to the 00:00 epochs and mostly commonly occur
    %  at 59:44. Therefore, if these are present, we want to use the new info
    %  and discard the ephemeris following it at 00:00.
    newuploads = find( mod( sat_ephem(:,17), 3600 ) ~= 0 );

    if ~isempty(newuploads) % If new uploads exist
        badrows = [];
        for zz = 1:length(newuploads)  % Loop over new uploads
            if (newuploads(zz) == length(sat_ephem(:,1))) % If new upload is last ephem entry
                continue;
            else
            % Check that new upload is within 4 minutes of next ephem entry
            %  In which case, remove the usual 00:00 entry and use the 
            %  new upload
            if ((sat_ephem(newuploads(zz)+1,17) - sat_ephem(newuploads(zz),17)) < 240)
                badrows = [badrows; kk(newuploads(zz) + 1)];    
            end
        end
    end % END of FOR loop
    % Remove unwanted ephemeris entries
    gps_ephem(badrows,:) = [];
    end % END of IF new uploads exist
    end
    end
end



end




