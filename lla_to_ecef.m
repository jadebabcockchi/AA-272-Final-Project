function ecef = lla_to_ecef(phi,lambda,h)
    phi = phi*(pi/180);
    lambda = lambda*(pi/180);
    Re = 6378.137*1000;         %[m] radius of the earth
    f = 1/298.257223563;        %flattening parameter of ellipsoid
    eSqrd = 2*f - f^2;          %square of eccentricity of elliosoid
    Ce = Re/sqrt(1-eSqrd*(sin(phi)^2)); %curvature of earth at some location
    
    x = (Ce + h)*cos(phi)*cos(lambda); %from lecture notes!
    y = (Ce + h)*cos(phi)*sin(lambda);
    z = ((Ce*(1 - eSqrd) + h)*sin(phi));
    ecef = [x y z];
end