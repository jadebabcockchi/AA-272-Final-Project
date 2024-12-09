function [az, el, range] = compute_az_el_range(user_ecef, sat_ecef)
    los_enu = compute_los_enu(user_ecef,sat_ecef); %unit vector pointing from user to satellite
    e = los_enu(1); n = los_enu(2); u = los_enu(3);
    az = atan2(e,n); 
    az = rad2deg(az);
    el = atan2(u,sqrt(e^2+n^2));
    el = rad2deg(el);
    range = norm(sat_ecef - user_ecef);
    %output in degrees
end