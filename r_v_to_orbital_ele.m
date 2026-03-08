function [a,e,i,RAAN,ap,theta] = r_v_to_orbital_ele(r0,v0,mu0, output_results)
    if nargin < 8
        output_results = false;
    end

    h = cross(r0,v0);
    n_hat = cross([0 0 1],h);
    e_vec = ((norm(v0)^2 - mu0/norm(r0)) * r0 - dot(r0,v0)*v0)/mu0;

    e = norm(e_vec);
    E = norm(v0)^2/2 - mu0/norm(r0);
    a = -mu0/(2*E);
    p = a*(1 - e^2);
    i = acosd(h(3)/norm(h));

    RAAN = acosd(n_hat(1)/norm(n_hat));
    if n_hat(2) < 0 % if y component of node is < 0
        RAAN = 360 - RAAN;
    end

    ap = acosd(dot(n_hat,e_vec) / (norm(n_hat)*e));
    if e_vec(3) < 0 % if z component < 0
        ap = 360 - ap;
    end

    theta = acosd(dot(e_vec,r0) / (e*norm(r0)));
    if dot(r0,v0) < 0
        theta = 360 - theta;
    end

    if output_results
        fprintf("a = %.4f\n",a);
        fprintf("p = %.4f\n",p);
        fprintf("i = %.4f deg\n",i);
        fprintf("Ω(longitude of ascending node) = %.4f deg\n",RAAN);
        fprintf("ω(arg. of periapses) = %.4f deg\n",ap);
        fprintf("θ_0(true anomaly) = %.4f deg\n",theta);
    end
end