%% MAE 569 Spring 2026 Final Project
% Contributers: Joeseph Presa, Lee Wallenfang, Ryan O'Connor, Devan McGaha

% This MATLAB script calls all other necessary scripts in order to solve the problem.
% Ensure that all necessary .m files are present in the directory when executing this script.

clear;clc;close all
% lower and upper bounds
% note: strange numbers in spot 1 are datenum format for start and end of
% launch window i.e. 2026-1-1, 2046-12-31

lb = [739983,... % launch time
    1, 1, 10, 60,... % transfer time, TU
    2000, 1500, 1000,... % flyby height, km
    0,0,0,0]; % path type

ub = [747652,... % launch time
    10, 10, 40, 475,... % transfer time, TU
    1e5, 1e5, 3e5,... % flyby height, km
    1,1,1,1]; % path type
 
% max expcected TU for hohmann transfers:
% e <-> m: 8.907 TU
% m-> j: 38.732 TU
% j-> n: 465.9 TU

% differential evolution options
opts.pop_size = 100;
opts.max_gen  = 100;
opts.F        = 0.8;
opts.CR       = 0.9;

fitness_fn = @(x) fitness(x);
[x_best, dv_best, x_best_list, dv_best_list] = differential_evolution(fitness_fn, lb, ub, opts);

vel_conversion_factor = 29.78594; % 1 AU/TU = 29.78594 km/s. dv is initially in AU/TU

% output evolution to excel doc
filename = "evodata.xlsx";
T = table(dv_best_list, x_best_list);
writetable(T, filename,'Sheet', 1);

% print basic results
format shortg
fprintf('Best solution: %.3f km/s total ΔV\n', dv_best*vel_conversion_factor);
fprintf('best state:\n');
fprintf("launch date: %s\n", datetime(x_best(1),'ConvertFrom','datenum'));
fprintf("earth to mars TOF: %f days\n", x_best(2) * 58.13);
fprintf("short/long way transfer? %i\n", x_best(9));
fprintf("mars flyby height: %f km\n", x_best(6));
fprintf("mars to earth TOF: %f days\n", x_best(3) * 58.13);
fprintf("short/long way transfer? %i\n", x_best(10));
fprintf("earth flyby height: %f km\n", x_best(7));
fprintf("earth to jupiter TOF: %f days\n", x_best(4) * 58.13);
fprintf("short/long way transfer? %i\n", x_best(11));
fprintf("jupiter flyby height: %f km\n", x_best(8));
fprintf("jupiter to neptune TOF: %f days\n", x_best(5) * 58.13);
fprintf("short/long way transfer? %i\n", x_best(12));

figure
hold on
plot(dv_best_list(:,1),dv_best_list(:,2))
hold off

% Polish with gradient descent
% options = optimset('MaxFunEvals', 50000, 'TolFun', 1e-6);
% [x_polished, dv_polished] = fminsearch(fitness_fn, x_best, options);
% 
% fprintf('true best: %.3f km/s total ΔV\n', dv_polished);
% fprintf('best state:');
% disp(x_polished);
