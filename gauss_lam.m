function [V1,V2] = gauss_lam(R1, R2, TOF, path_type)
    
    %   inputs:
    %    R1    - initial position vector [3x1], DU
    %    R2    - final position vector [3x1], DU
    %    TOF       - time of flight, TU
    %    path_type - 0 is short/propgrade, 1 is long/retrograde
    %
    %   outputs:
    %    v1 - vel at r1, DU/TU
    %    v2 - vel at r2, DU/TU
    
    mu = 1.0;
    tol = 1e-8;
    N_max = 500;
    
    r1 = norm(R1);
    r2 = norm(R2);
    theta = acos(dot(R1, R2) / (r1 * r2));
    if path_type == 1 % if path type is *long*
        theta = 2*pi - theta;
    end
    A = (sqrt(r1*r2)*sin(theta)) / sqrt(1 - cos(theta));

    
    % initial output
    % fprintf('  solving for %s path\n', path_type);
    % fprintf('  |r1| = %.6f DU   |r2| = %.6f DU   TOF = %.6f TU\n\n', r1, r2, TOF);

    % initial table row
    % fprintf('  %4s  %12s  %12s  %12s  %14s  %12s\n', ...
    %         'Iter', 'z', 'y', 'x', 't(z)', 'dt/dz');
    % fprintf('%s\n', repmat('=',1,72));

    z = 0.0;
    iter = 0;
    while true
        iter = iter + 1;

        % karl stumpff goated
        [C, S] = stumpff(z);

        y = r1 + r2 - A*(1 - z*S)/sqrt(C);

        x = sqrt(y/C);
        t_z = (x^3*S + A*sqrt(y)) / sqrt(mu);
        dtdz = get_dtdz(x, y, z, C, S, A, mu);

        dz = (TOF - t_z) / dtdz; % newton step

        
        % fprintf('  %4d  %12.6f  %12.6f  %12.6f  %14.8f  %12.6f\n', ...
                % iter, z, y, x, t_z, dtdz); % output results

        % break if converged
        if abs(TOF - t_z) < tol
            % fprintf('  Converged in %d iterations.  residual = %.3e TU\n\n', ...
                    % iter, abs(TOF - t_z));
            break
        end

        z_new = z + dz;
        z = z_new;

        % break and warn if max iter reached
        if iter >= N_max
            % warning('gauss_lambert: maximum iterations reached. Result may be inaccurate.');
            break
        end
    end

    % cleaning up results into v1, v2
    f = 1 - y/r1;
    g = A * sqrt(y/mu);
    fdot = ( -sqrt(mu)*x*(1 - z*S) ) / (r1*r2);
    gdot = 1 - y/r2;

    V1 = (R2 - f*R1) / g;
    V2 = (gdot*R2 - R1) / g;

    % lagrange identity sanity check
    lagrange_err = abs(f*gdot - fdot*g - 1);
    if lagrange_err > 1e-6
        % warning('Lagrange identity error = %.2f (check inputs)', lagrange_err);
    end

end


function [C, S] = stumpff(z)
    if z > 1e-6 % elliptic
        sqz = sqrt(z);
        C = (1 - cos(sqz))  / z;
        S = (sqz - sin(sqz)) / (z*sqz);
    elseif z < -1e-6 % hyperbolic
        sqz = sqrt(-z);
        C = (1 - cosh(sqz)) / z;
        S = (sinh(sqz) - sqz) / ((-z)*sqz);
    else % near-parabolic
        C = 0.5;
        S = 1/6;
        term = 1;
        for k = 1:20
            term = term * (-z) / ((2*k+2)*(2*k+3));
            S = S + term;
            % term_c = (-z)^k / factorial(2*k+2);
            C = C + (-1)^k * z^k / factorial(2*k+2);
            if abs(term) < 1e-15, break; end
        end
    end
end


function dtdz = get_dtdz(x, y, z, C, S, A, mu)
    % packed into its own function because what a hellish formula
    if abs(z) > 1e-6
        dCdz =  (1 - z*S - 2*C) / (2*z);
        dSdz =  (C - 3*S) / (2*z);
    else % there are really never enough edge cases
        dCdz = -1/24;
        dSdz = -1/120;
    end

    dtdz = (x^3*(dSdz - 3*S*dCdz/(2*C)) + ...
            A/8*(3*S*sqrt(y)/C + A/x)) / sqrt(mu);
end
