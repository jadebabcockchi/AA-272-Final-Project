function readRinex_obs_2024_GAL(obsfilename)

% pre-allocating an array 86400 rows (1 day of 1Hz data) 
% pre-allocating an array 7 columns (tow, gpsweek, prn, pseudorange,
% carrier phase, doppler, C/No)
obs = NaN(864000,7);

%% Header

%==========================================================================
headerend = [];
if (exist(obsfilename,'file') == 2)
    fid = fopen(obsfilename,'r');
else
    error(sprintf('Unable to find observation file: %s',obsfilename), 'ERROR!');
end

while (isempty(headerend) == 1)
   tline = fgetl(fid); 
   headerend = findstr(tline, 'END OF HEADER'); 
end

disp('end of header')
%==========================================================================

%% Body

j = 1;
while 1
    % Load next line in ephemeris file
    current_line = fgetl(fid);
    %disp(current_line)
    
    % If the next line is not a character then the end of the file has been
    %   reached and the while loop is exited
    if ~ischar(current_line)
        break; 
    end

    if current_line(1) == '>'

        yr = str2num(current_line(3:6));
             
        % Get the time for this data epoch.
        current_time = [ yr ; str2num(current_line(8:9)) ; ...
            str2num(current_line(11:12)) ; str2num(current_line(14:15)) ; ...
            str2num(current_line(17:18)) ; str2num(current_line(20:29)) ]';
        [gpswk, gpssec] = cal2gps2024(datetime(current_time));
    end
        
    if current_line(1) == 'E'
        if (~isempty(str2num(current_line(70:81)))) %(~isempty(str2num(current_line(21:34)))) % || 
            prn = str2num(current_line(2:3));
            pseudorange = str2num(current_line(70:81));
            cNo = str2num(current_line(124:129));

            if isempty(str2num(current_line(105:113)))
                doppler = NaN;
            else
                doppler = str2num(current_line(105:113));
            end
            if isempty(str2num(current_line(86:97)))
                carrier = NaN;
            else
                carrier = str2num(current_line(86:97));
            end
            obs(j,:) = [gpswk, gpssec, prn, pseudorange, carrier, doppler, cNo];
            
        else
            obs(j,1) = gpswk;
            obs(j,2) = gpssec;
        end
        j = j + 1;
    end

    % if mod(gpssec,1000)
    %     disp(gpssec)
    % end
end

matFilename = [obsfilename(1:end-4), '_GAL.mat'];
save(matFilename, 'obs');

end
