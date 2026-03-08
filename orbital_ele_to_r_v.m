function [R,V] = orbital_ele_to_r_v(a,e,i,Omega,ap,theta,mu, output_results)
    if nargin < 8
        output_results = false;
    end

    p = a * (1-e^2);
    h = sqrt(p * mu);

    rx = h^2 / mu * (1/(1 + e*cosd(theta))) * [cosd(theta); sind(theta); 0];
    vx = mu/h * [-sind(theta); (e + cosd(theta)); 0];

    rotation_mat = [cosd(ap), sind(ap),0;-sind(ap),cosd(ap),0;0,0,1]*...
    [1,0,0;0,cosd(i),sind(i);0,-sind(i),cosd(i)]*...
    [cosd(Omega), sind(Omega),0;-sind(Omega),cosd(Omega),0;0,0,1];

    R = rotation_mat\rx;
    V = rotation_mat\vx;

    if output_results
        disp(R)
        disp(V)
    end
end