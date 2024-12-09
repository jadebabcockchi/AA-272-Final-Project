function [R_exp, rGPS, bsv, relsv, all_el, all_az] = expectedRange_GAL(c, rRX, wn, tow, t_input, gal_ephem)
    %step 1: compute gps satellite postion in ECEF based on broadcast ephem.
    
    guess = rRX;
    for i = 1:36
        [~,pos,bsv(i),relsv(i)] = broadcast_eph2pos_modified(gal_ephem,t_input,i);
        rGPS(i,1) = pos(1); %x pos of satellite in ECEF
        rGPS(i,2) = pos(2); %y pos of satellite in ECEF
        rGPS(i,3) = pos(3); %z pos of satellite in ECEF
    end

    user_ecef = guess;
    for i = 1:36
        sat_ecef = rGPS(i,:);
        [az, el, range] = compute_az_el_range(user_ecef, sat_ecef);
        if i == 1
            all_az = az;
            all_el = el;
            all_range = range;
        else
            all_az = [all_az; az];
            all_el = [all_el; el];
            all_range = [all_range; range];
        end
    end

    %step 2: find range between RX and SAT with vector subtraction 

    Ro = rGPS - rRX; %GEOMETRIC RANGE
    for i = 1:length(Ro)
        R_old(i,1) = norm(Ro(i,:));
    end
    R_new = zeros(length(R_old),1);
    count = 0;

    diff = R_new - R_old;
    noNaN = diff(~isnan(diff));

    %step 3: compute time of transmission

    while (norm(noNaN) > 1e-5)
        clear noNaN
        if count > 0
            R_old = R_new;
        end
        Tr = tow.*ones(36,1);
        Tt = Tr - R_old./c;

        %step 4: recompute satellite position at Tt in ECEF based on signal travel time

        wn = wn.*ones(36,1);
        t_input = [wn Tt];
        for i = 1:36
            [~,pos,bsv(i),relsv(i)] = broadcast_eph2pos_modified(gal_ephem,t_input(i,:),i);
            rGPS_adj1(i,1) = pos(1);
            rGPS_adj1(i,2) = pos(2);
            rGPS_adj1(i,3) = pos(3);
        end

        %step 5: rotate satellite positions to ECEF at Tr

        wE  = 7.2921151467e-5; % WGS-84 value, rad/s
        for i = 1:length(rGPS_adj1) 
            phi = wE*(Tr-Tt);
            C_rot = [cos(phi(i)) sin(phi(i)) 0;...
                    -sin(phi(i)) cos(phi(i)) 0;...
                        0        0     1];
            if i == 1
                rGPS_adj2 = (C_rot*rGPS_adj1(i,:)')';
            else
                rGPS_adj2 = [rGPS_adj2; (C_rot*rGPS_adj1(i,:)')'];
            end
        end

        %step 6: compute new geometric range using this position for rGPS

        Rn = rGPS_adj2 - rRX;
        for i = 1:length(Rn)
            R_new(i,1) = norm(Rn(i,:));
        end
        count = count + 1;
        diff = R_new - R_old;
        noNaN = diff(~isnan(diff));
    end
    %repeat steps 3-6
    R_exp = R_new;
    
 end
