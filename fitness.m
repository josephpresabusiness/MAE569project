function dv_total = fitness (x)
    % outputs total delta v given state vector defining trajectory.
    % note: TOFs in tu, flyby heights in km
    % note: scalar r,v in lowercase, vector R,V in uppercase.

    %% unpacking state vector
    launch_time = x(1); % datetime after J2000
    TOFs = x(2:5); % e>m, m>e, e>j, j>n respectively. in TU
    flyby_heights = x(6:8); %m, e, j flyby heights respectively
    gauss_maneuver_directions = x(9:12); % 'short'/'long'
    % note for me: 2000 km orbit at neptune, 500 km orbit at earth.

    % ASSUMPTION MADE: time spent in flyby is negligible
    t1 = launch_time;
    t2 = launch_date + TOFs(1);
    t3 = t2 + TOFs(2);
    t4 = t3 + TOFs(3);
    t5 = t4 + TOFs(4);

    %% planet states + data
    [R_e0, V_e0] = r_v_after_j2000('Earth',   t1);
    [R_m,  V_m]  = r_v_after_j2000('Mars',    t2);
    [R_e1, V_e1] = r_v_after_j2000('Earth',   t3);
    [R_j,  V_j]  = r_v_after_j2000('Jupiter', t4);
    [R_n,  V_n]  = r_v_after_j2000('Neptune', t5);

    dv_total = 0;

    % TODO set actual data. r in km, mu in km^3/s^2
    rad_earth = 1;
    mu_earth = 1;
    rad_mars = 1;
    mu_mars = 1;
    rad_jup = 1;
    mu_jup = 1;
    rad_neptune = 1;
    mu_neptune = 1;

    %% leg 1: earth orbit to mars 
    
    v0 = sqrt(mu_earth/(rad_earth + 500)); % earth orbit velocity
    % TODO convert from km/s to AU/TU
    [V1, V2] = gauss_lam(R_e0, R_m, TOFs(1), gauss_maneuver_directions(1));
    % DV 1: leaving earth orbit
    dv_total = dv_total + norm(V1 - V_e0) - v0; % subtracting v0 to assume optimal launch angle

    %% leg 2: mars flyby
    % TODO write flyby func, make it assume ballistic trajectory. ignore
    % the horrifying mishmash of units in the mockup so far
    V_marsfb_out = flyby_exit_vel(V2, V_m, flyby_heights(1), rad_mars, mu_mars);

    %% leg 3: mars back to earth
    [V3, V4] = gauss_lam(R_m, R_e1, TOFs(2), gauss_maneuver_directions(2));
    % DV 2: DSM after mars flyby
    dv_total = dv_total + norm(V3 - V_marsfb_out);

    %% leg 4: earth flyby
    V_earthfb_out = flyby_exit_vel(V4, V_e1, flyby_heights(2), rad_earth, mu_earth);

    %% leg 5: earth to jupiter
    [V5, V6] = gauss_lam(R_e1, R_j, TOFs(3), gauss_maneuver_directions(3));
    % DV 3: DSM after earth flyby
    dv_total = dv_total + norm(V5 - V_earthfb_out);

    %% leg 6: jupiter flyby
    V_jupfb_out = flyby_exit_vel(V6, V_j, flyby_heights(3), rad_jup, mu_jup);

    %% leg 7: jupiter to neptune
    [V7, V8] = gauss_lam(R_j, R_n, TOFs(4), gauss_maneuver_directions(4));
    % DV 4: DSM after jupiter flyby
    dv_total = dv_total + norm(V7 - V_jupfb_out);

    %% leg 8: entry into neptune orbit
    v9 = sqrt(mu_neptune/(rad_neptune + 2000)); % neptune orbit velocity
    % DV 5: entering neptune orbit
    dv_total = dv_total + norm(V8 - V_n) - v9; % subtracting v9 to assume optimal entry angle
   
end