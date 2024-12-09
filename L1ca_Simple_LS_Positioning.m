% L1C/A Least Squares Positioning

clear; close all; clc;

f = 1;

if f == 1

    
    navfilename = 'L1L2mb04_ubxlog_3_all.nav'; %5845.3 5380 1257380
    obsfilename = 'L1L2mb04_ubxlog_3_all.obs';
    truth1 = ecef_to_lla(2176877.7800,621906.9814,5942863.6590);
    truth1ecef = [2176877.7800,621906.9814,5942863.6590];
    sample1 = [69.2730616018062	15.9439788454853	49.6096664518118]; % 7h 10m
    truth_lla = truth1;
    truth_ecef = truth1ecef;

elseif f == 2

    navfilename = 'L1L2mb05_ubxlog_4_all.nav'; %5379.8 5854 .... 6051400
    obsfilename = 'L1L2mb05_ubxlog_4_all.obs';
    truth2 = ecef_to_lla(2176681.9919,622095.7950,5942910.4073);
    truth2ecef = [2176681.9919,622095.7950,5942910.4073];
    sample2 = [69.2743635777591	15.9499449162027	47.3957500476390];
    truth_lla = truth2;
    truth_ecef = truth2ecef;

elseif f == 3

    navfilename = 'L1L5mb01_ubxlog_13_all.nav'; %5124.9 5125 .... 83000
    obsfilename = 'L1L5mb01_ubxlog_13_all.obs';
    truth3 = ecef_to_lla(2176776.6563,622527.8168,5942839.3876);
    truth3ecef = [2176776.6563,622527.8168,5942839.3876];
    sample3 = [69.2723750539215	15.9597245925882	54.8821056010202];
    truth_lla = truth3;
    truth_ecef = truth3ecef;

elseif f == 4
     
    navfilename = 'L1L5mb02_ubxlog_19_all.nav'; %ignore
    obsfilename = 'L1L5mb02_ubxlog_19_all.obs';
    truth4 = ecef_to_lla(2176289.7709,622770.3445,5942992.6023);
    truth4ecef = [2176289.7709,622770.3445,5942992.6023];
    sample4 = [69.2761797607958	15.9690816781175	37.7820149911568]; %9h 35m
    truth_lla = truth4;
    truth_ecef = truth4ecef;

else 

    navfilename = 'L1L5mb03_ubxlog_8_all.nav'; %5239.9 5240 .... 8879500
    obsfilename = 'L1L5mb03_ubxlog_8_all.obs';
    truth5 = ecef_to_lla(2176962.5171,622509.1930,5942777.2227);
    truth5ecef = [2176962.5171,622509.1930,5942777.2227];
    sample5 = [69.2707184490091	15.9579984489903	60.9691980155185];
    truth_lla = truth5;
    truth_ecef = truth5ecef;

end

%readRinex_obs_2024(obsfilename)



%% Jammertest in Bliek, Norway - Oslo, Norway Guess

%iterate your solution 5 times
%print out the new position and the provided true position

obs_data = [obsfilename(1:end-4), '.mat'];
load(obs_data);
[gps_ephem, ionoGPS] = read_clean_GPSbroadcast(navfilename,1);

zd = 2.375;
c   = 2.99792458e8;                                     %[m/s] SI speed of light
Jammertest_lla = [69.27 15.96 52]; 
lat = Jammertest_lla(1); lon = Jammertest_lla(2);
Oslo_lla = [59.9257 10.7466 0]; % make sure to guess something not too far away.... obs and all_el won't line up...
Oslo_ecef = lla2ecef(Oslo_lla);

% wn = 2331;
% tow = 210000; %201660;
% t_input = [wn tow];

hr = 7;
min = 0; %7h0m %7h10m
tow_start = min*60 + hr*3600 + 86400*2; %198859; 
hr = 7;
min = 0; %7h10m %8h5m
tow_end = min*60 + hr*3600 + 86400*2;
wn = 2331;

tows = tow_start:tow_end;

t_vec = [wn*ones(length(tows),1) tows'];
latlonalt_mat = NaN(length(tows),3);
latlonalt_ecef = NaN(length(tows),3);
clk_bias_mat = NaN(length(tows),1);
for t = 1:length(tows)
    clear A
    clear dy
    clear PIF
    clear R
    clear BSV
    clear RELSV
    clear IONO
    clear TROPO
    clear RGPS
    clear PSEUDO
    clear channels

    t_input = t_vec(t,:);
    wn = t_input(1,1);
    tow = t_input(1,2);


iter = 1;
guess_RX_loc = Oslo_ecef;
clock_guess = 0;

if t == 108
    disp('ok')
end

while iter < 6
    if iter == 1 && t == 3
        disp(iter)
    end
    %find expected range
    %clock error
    %elevation angles
    %azimuth angles
    [R_exp, rGPS, bsv, relsv, all_el, all_az] = expectedRange(c, guess_RX_loc, wn, tow, t_input, gps_ephem);
    epoch = tow;
    ionoklobuchar = klobuchar(ionoGPS, epoch, lat, lon, all_az, all_el);

    clear small_obs
    clear dy
    clear A
    clear PIF
    clear R
    clear BSV
    clear RELSV
    clear IONO
    clear TROPO
    clear RGPS
    clear PSEUDO
    clear channels

    % Logical indexing to find rows where column 2 is tow
    tow_rows = obs(:, 2) == tow;
    obs_prns = obs(tow_rows, 3);
    small_obs = obs(tow_rows,:);

    n = 1;
    for i = 2:32 %PRN 1 is out of commission right now
        if all_el(i) > 10
           index(n) = i;
           n = n + 1;
        end
    end
    shared = intersect(index, obs_prns);
    channels = length(shared);

    % for i = 1:length(obs_prns)
    %     small_obs = obs(tow_rows,:);
    % end

    if channels < 4
        guess_RX_loc = [NaN NaN NaN];
        clock_guess = NaN;
        fprintf("Jammertest ECEF Position: [%.2f, %.2f, %.2f] [m]\n",guess_RX_loc(1),guess_RX_loc(2),guess_RX_loc(3));
        break
    end

    for i = 1:channels
        prn = shared(i);
        R(i,1) = R_exp(prn);
        BSV(i,1) = bsv(prn);
        RELSV(i,1) = relsv(prn);
        IONO(i,1) = ionoklobuchar(prn);
        TROPO(i,1) = tropomodel(zd, all_el(prn));
        RGPS(i,:) = rGPS(prn,:);
        EL(i,1) = all_el(prn);
        AZ(i,1) = all_az(prn);
        row_pseudo = (small_obs(:, 2) == tow) & (small_obs(:, 3) == prn);
        PSEUDO(i,1)= small_obs(row_pseudo, 4);
    end
    
    A = zeros(channels,4); %geometry matrix
    for i = 1:channels
        A(i,1) = -(RGPS(i,1) - guess_RX_loc(1))/R(i);
        A(i,2) = -(RGPS(i,2) - guess_RX_loc(2))/R(i);
        A(i,3) = -(RGPS(i,3) - guess_RX_loc(3))/R(i);
        A(i,4) = 1;
    end
    
    dy = PSEUDO - (R - BSV - RELSV + IONO + TROPO); %measured values matrix
    dx = inv(A'*A)*A'*dy; %errors matrix
   
    guess_RX_loc = guess_RX_loc + dx(1:3)'; %change our guess and add the correction
    clock_guess = clock_guess + dx(4);
    
    %fprintf("Iteration: %d\n",iter);
    %fprintf("Position Errors: [%.2f, %.2f, %.2f] [m]\n",dx(1),dx(2),dx(3));
    %fprintf("Jammertest Adjusted ECEF Position: [%.2f, %.2f, %.2f] [m]\n",guess_RX_loc(1),guess_RX_loc(2),guess_RX_loc(3));
    %fprintf("NIST ECEF Position: [%.2f, %.2f, %.2f] [m]\n\n",NIST_adj(1),NIST_adj(2),NIST_adj(3));
    
    iter = iter + 1;
end

latlonalt = ecef_to_lla(guess_RX_loc(1),guess_RX_loc(2),guess_RX_loc(3))
clk_bias_mat(t,1) = clock_guess;
latlonalt_mat(t,:) = latlonalt;
latlonalt_ecef(t,:) = [guess_RX_loc(1),guess_RX_loc(2),guess_RX_loc(3)];

%check = ecef_to_lla(2176877.7800, 621906.9814, 5942863.6590)
%check = ecef_to_lla(2176681.9919, 622095.7950, 5942910.4073)

end

clk_sec = clk_bias_mat./c;
clk_sec = clk_sec - clk_sec(1);
clk_ns = clk_sec*10^9;

% wgs84 = wgs84Ellipsoid('meter');
% [xEast,yNorth,zUp] = ecef2enu(x,y,z,lat0,lon0,h0,wgs84)

%% Figures

% C_ecef_to_enu = ecef_to_enu(truth_lla(1),truth_lla(2));
% truth_enu = C_ecef_to_enu*truth_ecef';
% 
% for i = 1:length(latlonalt_mat)
%     pos_enu(:,i) = C_ecef_to_enu*latlonalt_ecef(i,:)';
% end
% 
% tow_values = tow_start:tow_end;
% hr_tod = (tow_values/86400 - 2)*24; % time of day [hr]
% min_toh = (hr_tod-7)*60; % time of hour ... in this case it's 07:00 UTC [min]
% ymin = -40;
% ymax = 40;
% 
% % plot postion residuals ENU
% figure()
% plot(min_toh, pos_enu(1,:)-truth_enu(1));
% hold on
% plot(min_toh, pos_enu(2,:)-truth_enu(2));
% plot(min_toh, latlonalt_mat(:,3)-truth_lla(3));
% xlim([0,10])
% ylim([-15,15])
% legend("E-error", "N-error", "U-error");
% xlabel("Time since 07:00 UTC 9/10/24 [min]");
% ylabel("Residual [m]");
% title("L1C/A RX5 Position Residuals vs. Time No RFI");
% grid on
% hold off
% 
% % figure()
% % plot(min_toh, pos_enu(1,:)-truth_enu(1));
% % hold on
% % plot(min_toh, pos_enu(2,:)-truth_enu(2));
% % plot(min_toh, latlonalt_mat(:,3)-truth_lla(3));
% % xlim([10,65])
% % ylim([ymin,ymax])
% % xShade1 = [15, 20, 20, 15];
% % xShade2 = [20, 25, 25, 20];
% % xShade3 = [35, 40, 40, 35];
% % xShade4 = [50, 55, 55, 50];
% % xShade5 = [55, 60, 60, 55];
% % yShade = [ymin, ymin, ymax, ymax];
% % fill(xShade1, yShade, 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
% % fill(xShade2, yShade, 'y', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
% % fill(xShade3, yShade, 'y', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
% % fill(xShade4, yShade, 'r', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
% % fill(xShade5, yShade, 'y', 'FaceAlpha', 0.05, 'EdgeColor', 'none');
% % legend("E-error", "N-error", "U-error");
% % xlabel("Time since 07:00 UTC 9/10/24 [min]");
% % ylabel("Residual [m]");
% % title("L1C/A RX5 Position Residuals vs. Time (Jamming + Meaconing)");
% % grid on
% % hold off
% % 
% % % splot clock bias
% % figure()
% % plot(min_toh,clk_sec-clk_sec(1));
% % hold on
% % grid on
% % xlabel("Time since 07:00 UTC 9/10/24 [min]");
% % ylabel("Clock Bias [s]")
% % title("Clock Bias vs. Time (Jamming + Meaconing)")
% % grid on
% % hold off



