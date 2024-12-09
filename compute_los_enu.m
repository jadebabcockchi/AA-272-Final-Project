function los_enu = compute_los_enu(user_ecef,sat_ecef)
    user_lla = ecef_to_lla(user_ecef(1),user_ecef(2),user_ecef(3)); %[m] LLA
    lat = user_lla(1);
    lon = user_lla(2);
    C_ecef_to_enu = ecef_to_enu(lat,lon);       %matrix to transform los_ecef to los_enu
    los_ecef = sat_ecef - user_ecef;            %find los vector in ecef frame
    [~,unit_los_ecef] = unitvec(los_ecef,2);    %normalize los_ecef
    los_enu = C_ecef_to_enu*unit_los_ecef';     %los_enu as a unit vector
end