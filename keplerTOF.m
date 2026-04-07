function [r_vec, v_vec] = keplerTOF(r0_vec,v0_vec,mu,t)
    r0 = norm(r0_vec);
    % v0 = norm(v0_vec);
    rdv = dot(r0_vec,v0_vec);

    [a,~,~,~,~,~] = r_v_to_orbital_ele(r0_vec,v0_vec,mu); % i am so clever
    sqra = sqrt(a);
    sqrma = sqrt(mu*a);
    sqrm = sqrt(mu);

    mu_tn = @(xn) a*( xn - sqra*sin(xn/sqra)) + (rdv*a*( 1 - cos(xn/sqra) )/sqrm) + ...
    (r0*sqra*sin(xn/sqra));
    tn = @(xn) (mu_tn(xn) / sqrm) - t;
    x = fzero(tn,abs(sqrma*t));

    %calculate f and g, get r
    f = 1 - a*(1 - cos(x/sqra))/r0;
    g = (a^2 / sqrma) * ( rdv*(1 - cos(x/sqra))/sqrma + r0*(sin(x/sqra))/a );
    r_vec = f*r0_vec + g*v0_vec;
    r = norm(r_vec);

    %calculate fdot, gdot, get v
    fdot = -sqrma*sin(x/sqra) / (r0*r);
    gdot = 1 - (a/r) + (a/r)*cos(x/sqra);
    v_vec = fdot*r0_vec + gdot*v0_vec;

    %sanity check
    % fprintf("%.4f should equal 1.\n",f*gdot - g*fdot)
end
