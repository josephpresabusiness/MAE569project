close all;
time = datetime(2026,12,25,11,12,0); % actual date

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
% ax.XColor = [1 1 1];       % white axis lines/labels
% ax.YColor = [1 1 1];
% ax.ZColor = [1 1 1];
hold on; axis equal; % grid on;
set(gcf, 'Color', 'k');

%% stars for fun
% size=40; % 40 TU spread of stars
% num_stars = 200;
% x=size*rand(1,num_stars); y=size*rand(1,num_stars); z=size*rand(1,num_stars);
% scatter3(2*(x-0.5*size), 2*(y-0.5*size), 2*(z-0.5*size), 'MarkerFaceColor', 'w', 'LineWidth', 0.1);

%% sun

scatter3([0],[0],[0], 'MarkerFaceColor', 'y', 'LineWidth', 1);

%% graph planet orbits
selected_planets = ["Earth", "Mars", "Jupiter", "Neptune", ...
    "Mercury", "Venus", "Saturn", "Uranus", "Pluto"];
planet_colors = ["#11BF59", "#BF3411", "#BF8B11", "#1173BF", ...
    "#F4740C", "#C9BE46", "#BCB993", "#22D9E2", "#B3BFD6"];
idx = 1;
num_points = 300;

for planet = selected_planets
    color = planet_colors(idx);
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

[R_me, ~] = r_v_after_j2000('Mercury',   time);
[R_ve,  ~]  = r_v_after_j2000('Venus',    time);
[R_ea, ~] = r_v_after_j2000('Earth',   time);
[R_ma, ~] = r_v_after_j2000('Mars',   time);
[R_ju,  ~]  = r_v_after_j2000('Jupiter', time);
[R_sa, ~] = r_v_after_j2000('Saturn',   time);
[R_ur, ~] = r_v_after_j2000('Uranus',   time);
[R_ne,  ~]  = r_v_after_j2000('Neptune', time);
[R_pl, ~] = r_v_after_j2000('Pluto',   time);
planet_positions = [R_ea, R_ma, R_ju, R_ne, R_me, R_ve, R_sa, R_ur, R_pl];

for i=1:9
    scatter3(planet_positions(1,i), planet_positions(2,i), planet_positions(3,i), ...
        120, 'filled', ...
        'MarkerFaceColor', planet_colors(i), 'DisplayName', selected_planets(i), 'Clipping', 'off');
    % Optional: label them
    % text(planet_positions(1,i), planet_positions(2,i), planet_positions(3,i), ...
    %     "  " + selected_planets(i), 'Color', 'w', 'FontSize', 9);
end

gca.Clipping = 'off';
