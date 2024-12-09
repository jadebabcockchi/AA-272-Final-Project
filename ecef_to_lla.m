function lla = ecef_to_lla(x,y,z)
    rho = sqrt(x^2+y^2);
    r = sqrt(x^2+y^2+z^2);
    oldGuess = 0;
    newGuess = asin(z/r);
    Re = 6378.137*1000;         %[m] radius of the earth
    f = 1/298.257223563;        %flattening parameter of ellipsoid
    eSqrd = 2*f - f^2;          %square of eccentricity of elliosoid

    %find the lattitude angle with a tolerance
    while (abs(oldGuess-newGuess) > 1e-8)
        oldGuess = newGuess;
        Ce = Re/sqrt(1-eSqrd*(sin(oldGuess)^2)); %curvature of earth at a location
        newGuess = atan((z + Ce*eSqrd*sin(oldGuess))/rho);
    end

    %define lat, lon, and height values for NIST
    phi = newGuess;
    lat = rad2deg(phi);
    lambda = atan2(y,x);
    lon = rad2deg(lambda);
    Ce = Re/sqrt(1-eSqrd*(sin(phi)^2));
    h = (rho/cos(phi))-Ce;
    lla = [lat lon h];
end