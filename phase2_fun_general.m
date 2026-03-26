function [V1,V2] = phase2_fun_general(TOF,planet1,planet2,date,path_type)

% Inputs:
% TOF: Time of Flight
% planet1: first planet
% planet2: second planet
% date: time since J2000 at first planet
% path_type: "short" or "long". designates prograde/retrograde trajectories
% respectively
%
% Outputs:
% V1: V_inf at planet 1 in km/s in Sun Inertial
% V2: V_inf at planet 2 in km/s in Sun Inertial
%
%
 global r1 r2 mu A
 location1 = prop(planet1,date);
 location2 = prop(planet2,date+TOF);
 R1 = location1;
 R2 = location2;
 AU = 149597870700/1000;
 TU = 58.13*24*60*60;
 VU = AU/TU;
 mu = 1;


 %% Magnitudes of R1 and R2
r1 = norm(R1);
r2 = norm(R2);

%% Theta 

c12 = cross(R1, R2);
if c12(3) <= 0
    theta = 2*pi - theta;
end

if path_type == "short"
    theta = acos(dot(R1,R2)/r1/r2); %prograde trajectory (short way solution
elseif path_type == "long"
    theta_short = acos(dot(R1,R2)/r1/r2); %retrograde trajectory (long way solution)
    theta=2*pi-theta_short;
else
    error("bad path type chosen. must be 'short'  or 'long'.")
end


%% Solve equation for z
A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));
%Check where F(z,t) changes sign, and use that 
%value of z as a starting point 
z = -100;
while F(z,t) < 0
    z = z + 0.1;
end
%Tolerence and number of iterations
eps= 1.e-8;
nmax = 5000;
%Iterate until z is found under the imposed tolerance:
ratio = 1;
n = 0;
while (abs(ratio) > eps) && (n <= nmax)
    n = n + 1;
    ratio = F(z,t)/dFdz(z);
    z = z - ratio;
end
%Maximum number of iterations exceeded
if n >= nmax
    fprintf('\n\n **Number of iterations exceeds %g \n\n ',nmax)
end

%% Lagrange coefficents
f = 1 - y(z)/r1;
g = A*sqrt(y(z)/mu);
g_dot = 1 - y(z)/r2;

%% V1 and V2
V1 = 1/g*(R2 - f*R1)*VU;
V2 = 1/g*(g_dot*R2 - R1)*VU;


%% Subfunctions used 
function dum = y(z)
    dum = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
end

function dum = F(z,t)
    dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
end

function dum = dFdz(z)
    if z == 0
        dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
    else
        dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
        + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...
        + A*sqrt(C(z)/y(z)));

    end
end

%...Stumpff functions:
function c = C(z)
    if z > 0
        c = (1 - cos(sqrt(z)))/z;
    elseif z < 0
        c = (cosh(sqrt(-z)) - 1)/(-z);
    else
        c = 1/2;
    end
end

function s = S(z)
    if z > 0
        s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
    elseif z < 0
        s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
    else
        s = 1/6;
    end
end
end