clc
format longG
planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];

target_date = datetime(2026,12,15,11,12,0); % actual date
% target_date = datetime(2023,12,1,5,11,0); % test date. anomalies of earth and mars should be 325.9 and 262.98

r_vec_out = zeros(num_planets,3);
v_vec_out = zeros(num_planets,3);
theta_out = zeros(num_planets,1);

for pl = 1:9
    planet = planets(pl);

    [R, V, theta] = r_v_after_j2000(planet, target_date);

    % assign data to output arrays
    r_vec_out(pl,:) = round(R,4);
    v_vec_out(pl,:) = round(V,4);
    theta_out(pl) = round(theta,4);
end

%% output data
output = table(planets',r_vec_out, v_vec_out, theta_out,'VariableNames',{'Planet', 'R(AU)', 'V(AU/TU)', 'True Anomaly (deg)'});
disp(output)
