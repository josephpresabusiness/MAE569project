function [dv,TOF] = gravity_assist(v1a,v1d,approach_alt,r_planet,v_planet,mu,i_planet1,i_planet2,m_planet1,R_planet1)
    % Inputs
    % v1a: approach velocity in sun frame
    % v1d: departure velocity in sun frame
    % approach alt: altitude of approach
    % r_planet: radius of planet
    % v_planet: velocity of planet in sun frame
    % mu: gravitational parameter
    % i_planet1: inclination of planet 1
    % i_planet2: inclinaton of planet 2

    % Outputs
    % dv: delta-v to go from v1a to v1d, including a burn at the turn point
    M_sun = 1.989*10^30;
    r_SOI = R_planet1*(m_planet1/M_sun)^2/5;
    rp = r_planet+approach_alt; % scalar
    v_inf_a_vec = v1a-v_planet; % vector
    v_inf_d_vec = v1d-v_planet; % vector
    v_inf_a = norm(v_inf_a_vec); % scalar
    v_inf_d = norm(v_inf_d_vec); % scalar
    epsilon_a = v_inf_a^2/2; % scalar
    epsilon_d = v_inf_d^2/2; % scalar
    a_a = -mu/(2*epsilon_a); % scalar
    a_d = -mu/(2*epsilon_d); % scalar
    b_a = sqrt((rp-a_a)^2-a_a^2); % scalar
    b_d = sqrt((rp-a_d)^2-a_d^2); % scalar
    e_a = sqrt(1+(b_a/a_a)^2); % scalar
    e_d = sqrt(1+(b_d/a_d)^2); % scalar
    h_a_vec = b_a*v_inf_a_vec; % vector
    h_d_vec = b_d*v_inf_d_vec; % vector
    h_a = norm(h_a_vec); % scalar
    h_d = norm(h_d_vec); % scalar
    vp_a = h_a/rp; % scalar
    vp_d = h_d/rp; % scalar
    delta_a = 2*asind(1/e_a); % degree
    delta_d = 2*asind(1/e_d); % degree
    delta_i = i_planet2-i_planet1; % degree
    dv = sqrt(vp_a^2 +vp_d^2 -2*vp_a*vp_d*cosd(delta_i)); % scalar
    n_a = sqrt(mu/((-a_a)^3));
    n_d = sqrt(mu/((-a_d)^3));
    p_a = h_a^2/mu;
    p_d = h_d^2/mu;
    theta_a = acosd(1/e_a*(p_a/r_SOI -1));
    theta_p = 0;
    theta_d = acosd(1/e_d*(p_d/r_SOI-1));
    F_a = acosh((e_a+cosd(theta_a))/(1+e_a*cosd(theta_a)));
    F_pa = acosh((e_a+cosd(theta_p))/(1+e_a*cosd(theta_p)));
    F_pd = acosh((e_d+cosd(theta_p))/(1+e_d*cosd(theta_p)));
    F_d = acosh((e_d+cosd(theta_d))/(1+e_d*cosd(theta_d)));
    TOF_a = ((e_a*sinh(F_a)-F_a)-(e_a*sinh(F_pa)-F_pa))/n_a;
    TOF_d = ((e_d*sinh(F_d)-F_d)-(e_d*sinh(F_pd)-F_pd))/n_d;
    TOF = TOF_a+TOF_d;
end
