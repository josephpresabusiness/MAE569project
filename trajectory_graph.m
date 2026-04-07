function [] = trajectory_graph(x)
    close all;
    which gca;
    
    %% unpacking state vector
    launch_time = x(1); % datetime after J2000
    TOFs = x(2:5); % e>m, m>e, e>j, j>n respectively. in TU
    % flyby_heights = x(6:8); %m, e, j flyby heights respectively (km)
    gauss_maneuver_directions = x(9:12); % 'short'/'long'
    
    %% flight times
    t1 = launch_time;
    t2 = launch_time + TOFs(1)*58.13;
    t3 = t2 + TOFs(2)*58.13;
    t4 = t3 + TOFs(3)*58.13;
    t5 = t4 + TOFs(4)*58.13;
    
    %% J2000 planet data
    planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];
    a_list = [0.387099 0.723332 1 1.523662 5.203363 9.53707 19.19126 30.06896 39.48169]; % AU
    e_list = [0.205631 0.006773 0.01671 0.093412 0.048393 0.054151 0.047168 0.008586 0.248808];
    i_list = [7.00487 3.39471 0.00005 1.85061 1.3053 2.48446 0.76986 1.76917 17.14175]; % deg
    Omega_list = [48.33167 76.68069 -11.26064 49.57854 100.55615 113.71504 74.22988 131.72169 110.30347]; % longitude of ascending node, deg
    ap_list = [29.12478 54.85229 114.20783 286.4623 -85.8023 -21.2831 96.73436 -86.75034 113.76329]; % arg of periapses, deg
    % theta_list = [174.7944 50.44675 -2.48284 19.41248 19.55053 -42.4876 142.2679 259.9087 14.86205]; % true anomaly, deg
    
    %% initialize figure
    f = figure; 
    ax = axes(f);
    f.Color = [0 0 0];
    ax.Color = [0 0 0];
    ax.XColor = [1 1 1];       % white axis lines/labels
    ax.YColor = [1 1 1];
    ax.ZColor = [1 1 1];
    hold on; axis equal; % grid on;
    set(gcf, 'Color', 'k');
    
    %% stars for fun
    % size=40; % 40 TU spread of stars
    % num_stars = 200;
    % x=size*rand(1,num_stars); y=size*rand(1,num_stars); z=size*rand(1,num_stars);
    % scatter3(2*(x-0.5*size), 2*(y-0.5*size), 2*(z-0.5*size), 'MarkerFaceColor', 'w', 'LineWidth', 0.1);
    
    
    %% graph planet orbits
    selected_planets = ["Earth", "Mars", "Jupiter", "Neptune"];
    orbit_colors = ["#11BF59", "#BF3411", "#BF8B11", "#1173BF"];
    idx = 1;
    num_points = 300;
    
    for planet = selected_planets
        color = orbit_colors(idx);
        r_traj = zeros(3,num_points);
        index = find(planets == planet); % index of planet in lists
        a = a_list(index); e = e_list(index); i = i_list(index);
        Om = Omega_list(index); ap = ap_list(index);
        theta_space = linspace(0,360,num_points);
    
        for k = 1:num_points
            [r_traj(:, k),~] = orbital_ele_to_r_v(a,e,i,Om,ap,theta_space(k),1,false);
        end
    
        plot3(r_traj(1,:), r_traj(2,:), r_traj(3,:), ...
            '--', 'Color', color, 'Clipping', 'off');
    
        idx = idx + 1;
    end
    
    view(30, 25);
    %% graph planet positions for each step
    
    [R_e0, ~] = r_v_after_j2000('Earth',   t1);
    [R_m,  ~]  = r_v_after_j2000('Mars',    t2);
    [R_e1, ~] = r_v_after_j2000('Earth',   t3);
    [R_j,  ~]  = r_v_after_j2000('Jupiter', t4);
    [R_n,  ~]  = r_v_after_j2000('Neptune', t5);
    planet_positions = [R_e0, R_m, R_e1, R_j, R_n];
    
    selected_planets = ["Earth", "Mars", "Earth(2nd)", "Jupiter", "Neptune"];
    planet_colors = ["#11BF59", "#BF3411", "#11BF59", "#BF8B11", "#1173BF"];
    for i=1:5
        scatter3(planet_positions(1,i), planet_positions(2,i), planet_positions(3,i), ...
            120, 'filled', ...
            'MarkerFaceColor', planet_colors(i), 'DisplayName', selected_planets(i), 'Clipping', 'off');
        % Optional: label them
        text(planet_positions(1,i), planet_positions(2,i), planet_positions(3,i), ...
            "  " + selected_planets(i), 'Color', 'w', 'FontSize', 9);
    end
    
    
    %% graph trajectories
    
    % earth to mars
    [V1, ~] = gauss_lam(R_e0, R_m, TOFs(1), gauss_maneuver_directions(1));
    % mars to earth
    [V3, ~] = gauss_lam(R_m, R_e1, TOFs(2), gauss_maneuver_directions(2));
    % earth to jupiter
    [V5, ~] = gauss_lam(R_e1, R_j, TOFs(3), gauss_maneuver_directions(3));
    % jupiter to neptune
    [V7, ~] = gauss_lam(R_j, R_n, TOFs(4), gauss_maneuver_directions(4));
    
    num_points = 200;
    traj_R = [R_e0, R_m, R_e1, R_j];
    traj_V = [V1,V3,V5,V7];
    for k = 1:4
        r_traj = zeros(3,num_points);
        TOF_space = linspace(0.01,TOFs(k),num_points);
        R = traj_R(:,k);
        V = traj_V(:,k);
    
        for i = 1:num_points
            [r_traj(:, i),~] = keplerTOF(R,V,1,TOF_space(i));
        end 
        disp(k);
        plot3(r_traj(1,:), r_traj(2,:), r_traj(3,:), ...
                '-', 'Color', 'w', 'Clipping', 'off');
    end
    
    gca.Clipping = 'off';
end
clc; clear all; close all;
% ignore this bit for now
x = [739983,... % launch time
    2, 5, 10, 60,... % transfer time, TU
    2000, 1500, 1000,... % flyby height, km
    1,0,1,0]; % path type
%% unpacking state vector
launch_time = x(1); % datetime after J2000
TOFs = x(2:5); % e>m, m>e, e>j, j>n respectively. in TU
flyby_heights = x(6:8); %m, e, j flyby heights respectively (km)
gauss_maneuver_directions = x(9:12); % 'short'/'long'

%% flight times
t1 = launch_time;
t2 = launch_time + TOFs(1)*58.13;
t3 = t2 + TOFs(2)*58.13;
t4 = t3 + TOFs(3)*58.13;
t5 = t4 + TOFs(4)*58.13;

%% J2000 planet data
planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];
a_list = [0.387099 0.723332 1 1.523662 5.203363 9.53707 19.19126 30.06896 39.48169]; % AU
e_list = [0.205631 0.006773 0.01671 0.093412 0.048393 0.054151 0.047168 0.008586 0.248808];
i_list = [7.00487 3.39471 0.00005 1.85061 1.3053 2.48446 0.76986 1.76917 17.14175]; % deg
Omega_list = [48.33167 76.68069 -11.26064 49.57854 100.55615 113.71504 74.22988 131.72169 110.30347]; % longitude of ascending node, deg
ap_list = [29.12478 54.85229 114.20783 286.4623 -85.8023 -21.2831 96.73436 -86.75034 113.76329]; % arg of periapses, deg
theta_list = [174.7944 50.44675 -2.48284 19.41248 19.55053 -42.4876 142.2679 259.9087 14.86205]; % true anomaly, deg

%% initialize figure
figure; hold on; axis equal; % grid on;
set(gca, 'Color', 'k'); % black background
set(gcf, 'Color', 'k');
ax = gca;
ax.XColor = 'w'; ax.YColor = 'w'; ax.ZColor = 'w';

%% stars for fun
% size=40; % 40 TU spread of stars
% num_stars = 200;
% x=size*rand(1,num_stars); y=size*rand(1,num_stars); z=size*rand(1,num_stars);
% scatter3(2*(x-0.5*size), 2*(y-0.5*size), 2*(z-0.5*size), 'MarkerFaceColor', 'w', 'LineWidth', 0.1);


%% graph planet orbits
selected_planets = ["Earth", "Mars", "Jupiter", "Neptune"];
orbit_colors = ["#11BF59", "#BF3411", "#BF8B11", "#1173BF"];
year_TU = 365.25636/58.13; % 1 year in TU
idx = 1;
num_points = 300;

for planet = selected_planets
    color = orbit_colors(idx);
    r_traj = zeros(3,num_points);
    index = find(planets == planet); % index of planet in lists
    a = a_list(index); e = e_list(index); i = i_list(index);
    Om = Omega_list(index); ap = ap_list(index);
    theta_space = linspace(0,360,num_points);

    for k = 1:num_points
        [r_traj(:, k),~] = orbital_ele_to_r_v(a,e,i,Om,ap,theta_space(k),1,false);
    end

    plot3(r_traj(1,:), r_traj(2,:), r_traj(3,:), ...
        '--', 'Color', color, 'Clipping', 'off');

    idx = idx + 1;
end

view(30, 25);
%% graph planet positions for each step

[R_e0, ~] = r_v_after_j2000('Earth',   t1);
[R_m,  ~]  = r_v_after_j2000('Mars',    t2);
[R_e1, ~] = r_v_after_j2000('Earth',   t3);
[R_j,  ~]  = r_v_after_j2000('Jupiter', t4);
[R_n,  ~]  = r_v_after_j2000('Neptune', t5);
planet_positions = [R_e0, R_m, R_e1, R_j, R_n];

selected_planets = ["Earth", "Mars", "Earth(2nd)", "Jupiter", "Neptune"];
planet_colors = ["#11BF59", "#BF3411", "#11BF59", "#BF8B11", "#1173BF"];
for i=1:5
    scatter3(planet_positions(1,i), planet_positions(2,i), planet_positions(3,i), ...
        120, 'filled', ...
        'MarkerFaceColor', planet_colors(i), 'DisplayName', selected_planets(i), 'Clipping', 'off');
    % Optional: label them
    text(planet_positions(1,i), planet_positions(2,i), planet_positions(3,i), ...
        "  " + selected_planets(i), 'Color', 'w', 'FontSize', 9);
end


%% graph trajectories

% earth to mars
[V1, ~] = gauss_lam(R_e0, R_m, TOFs(1), gauss_maneuver_directions(1));
% mars to earth
[V3, ~] = gauss_lam(R_m, R_e1, TOFs(2), gauss_maneuver_directions(2));
% earth to jupiter
[V5, ~] = gauss_lam(R_e1, R_j, TOFs(3), gauss_maneuver_directions(3));
% jupiter to neptune
[V7, ~] = gauss_lam(R_j, R_n, TOFs(4), gauss_maneuver_directions(4));

num_points = 200;
traj_R = [R_e0, R_m, R_e1, R_j];
traj_V = [V1,V3,V5,V7];
for k = 1:4
    r_traj = zeros(3,num_points);
    TOF_space = linspace(0.01,TOFs(k),num_points);
    R = traj_R(:,k);
    V = traj_V(:,k);

    for i = 1:num_points
        [r_traj(:, i),~] = keplerTOF(R,V,1,TOF_space(i));
    end 
    disp(k);
    plot3(r_traj(1,:), r_traj(2,:), r_traj(3,:), ...
            '-', 'Color', 'w', 'Clipping', 'off');
end

gca.Clipping = 'off';
