function C_ecef_to_enu = ecef_to_enu(lat_deg,lon_deg)
    phi = deg2rad(lat_deg);     %radians
    lambda = deg2rad(lon_deg);  %radians
    C = [-sin(lambda)           cos(lambda)             0; %from lecture notes!
        -sin(phi)*cos(lambda)   -sin(phi)*sin(lambda)   cos(phi)
        cos(phi)*cos(lambda)    cos(phi)*sin(lambda)    sin(phi)];
    C_ecef_to_enu = C;
end