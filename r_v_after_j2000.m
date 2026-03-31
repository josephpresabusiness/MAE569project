function [R2,V2,theta2] = r_v_after_j2000(planet,target_date)
    planets = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];
    a_list = [0.387099 0.723332 1 1.523662 5.203363 9.53707 19.19126 30.06896 39.48169]; % AU
    e_list = [0.205631 0.006773 0.01671 0.093412 0.048393 0.054151 0.047168 0.008586 0.248808];
    i_list = [7.00487 3.39471 0.00005 1.85061 1.3053 2.48446 0.76986 1.76917 17.14175]; % deg
    Omega_list = [48.33167 76.68069 -11.26064 49.57854 100.55615 113.71504 74.22988 131.72169 110.30347]; % longitude of ascending node, deg
    ap_list = [29.12478 54.85229 114.20783 286.4623 -85.8023 -21.2831 96.73436 -86.75034 113.76329]; % arg of periapses, deg
    theta_list = [174.7944 50.44675 -2.48284 19.41248 19.55053 -42.4876 142.2679 259.9087 14.86205]; % true anomaly, deg

    mu_sun = 1; % Au^3 / TU^2

    if isa(target_date, "double")
        target_date = datetime(target_date,'ConvertFrom','datenum');
    end
    J2000 = datetime(2000,1,1,11,58,0);
    dt = days(target_date - J2000)/58.13; % TU (1 TU in mu=1 space is 1/2pi years = 58.13 days)
    
    index = find(planets == planet); % index of planet in lists
    a = a_list(index);
    e = e_list(index);
    i = i_list(index);
    Om = Omega_list(index);
    ap = ap_list(index);
    theta = theta_list(index);

    % modulating dt by period of orbit
    period = (2*pi*sqrt(a^3/mu_sun)); % s
    dt_planet = mod(dt, period); 

    % convert to vector
    [R1, V1] = orbital_ele_to_r_v(a,e,i,Om,ap,theta,mu_sun);

    %get R2,V2 from kepler solver
    [R2, V2] = keplerTOF(R1,V1,mu_sun,dt_planet);

    % extract theta from output R2,V2 in the smartest way possible(lol)
    [~,~,~,~,~,theta2] = r_v_to_orbital_ele(R2,V2,mu_sun);
end
