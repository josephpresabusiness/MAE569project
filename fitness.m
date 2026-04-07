function dv_total = fitness (x)
    % outputs total delta v given state vector defining trajectory.
    % note: TOFs in tu, flyby heights in km
    % note: scalar r,v in lowercase, vector R,V in uppercase.

    %% convenient constants
    TU_s = 5022432; % 1 TU in seconds
    AU_km = 149597870.7; % 1 AU in km

    %% unpacking state vector
    launch_time = x(1); % datetime after J2000
    TOFs = x(2:5); % e>m, m>e, e>j, j>n respectively. in TU
    flyby_heights = x(6:8); %m, e, j flyby heights respectively (km)
    gauss_maneuver_directions = x(9:12); % 'short'/'long'
    % note for me: 2000 km orbit at neptune, 500 km orbit at earth.

    % ASSUMPTION MADE: time spent in flyby is negligible. this is a fairly
    % standard practice in preliminary mission design, so we feel fairly
    % safe in assuming it. 
    t1 = launch_time;
    t2 = launch_time + TOFs(1)*58.13;
    t3 = t2 + TOFs(2)*58.13;
    t4 = t3 + TOFs(3)*58.13;
    t5 = t4 + TOFs(4)*58.13;

    %% planet states + data

    [R_e0, V_e0] = r_v_after_j2000('Earth',   t1);
    [R_m,  V_m]  = r_v_after_j2000('Mars',    t2);
    [R_e1, V_e1] = r_v_after_j2000('Earth',   t3);
    [R_j,  V_j]  = r_v_after_j2000('Jupiter', t4);
    [R_n,  V_n]  = r_v_after_j2000('Neptune', t5);

    dv_total = 0;

    % radius and std gravitational parameter for relevant planets
    % TODO check these conversions bc i trust nothing

    rad_earth = 6378.1366 / AU_km;
    mu_earth = 398600 * (TU_s^2 / AU_km^3);
    rad_mars = 3396.19 / AU_km;
    mu_mars = 42828 * (TU_s^2 / AU_km^3);
    rad_jup = 71492 / AU_km;
    mu_jup = 126686000 * (TU_s^2 / AU_km^3);
    rad_neptune = 24764 / AU_km;
    mu_neptune = 6835100 * (TU_s^2 / AU_km^3); 

    %% calculating velocity vectors for transfer orbits
    % earth to mars
    [V1, V2] = gauss_lam(R_e0, R_m, TOFs(1), gauss_maneuver_directions(1));

    % mars to earth
    [V3, V4] = gauss_lam(R_m, R_e1, TOFs(2), gauss_maneuver_directions(2));

    % earth to jupiter
    [V5, V6] = gauss_lam(R_e1, R_j, TOFs(3), gauss_maneuver_directions(3));

    % jupiter to neptune
    [V7, V8] = gauss_lam(R_j, R_n, TOFs(4), gauss_maneuver_directions(4));
    

    %% deltaV from flybys
    % TODO add inclination params. zeroes currently

    % mars flyby
    dv_marsfb = gravity_assist(V2, V3, flyby_heights(1), rad_mars, V_m, mu_mars,1.85, 0.00005, 0.641691*10^24, R_m);

    % earth flyby
    dv_earthfb = gravity_assist(V4, V5, flyby_heights(2), rad_earth, V_e1, mu_earth, 0.00005, 1.304, 5.97217*10^24, R_e1);

    % jupiter flyby
    dv_jupfb = gravity_assist(V6, V7, flyby_heights(3), rad_jup, V_j, mu_jup, 1.304, 1.77, 1898.125*10^24, R_j);

    % add flyby dv accrued to total
    dv_total = dv_total + dv_marsfb + dv_earthfb + dv_jupfb;


    %% deltaV incurred from leaving/entering orbits

    % DV leaving earth orbit
    v0 = sqrt(mu_earth/(rad_earth + 500)) * (TU_s/AU_km); % earth orbit velocity
    dv_total = dv_total + norm(V1 - V_e0) - v0; % subtracting v0 to assume optimal launch angle
    
    % DV entering neptune orbit
    v9 = sqrt(mu_neptune/(rad_neptune + 2000)) * (TU_s/AU_km); % neptune orbit velocit
    dv_total = dv_total + norm(V8 - V_n) - v9 ; % subtracting v9 to assume optimal entry angle
   
end
