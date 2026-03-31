% lower and upper bounds
% note: strange numbers in spot 1 are datenum format for start and end of
% launch window i.e. 2026-1-1, 2046-12-31
lb = [739983,... % launch time
    1,  1,  10, 60,... % transfer time, TU
    200,   200,   1e5,... % flyby height, km
    0,0,0,0]; % path type
ub = [747652,... % launch time
    10, 10, 40, 475,... % transfer time, TU
    10000, 10000, 3e5,... % flyby height, km
    1,1,1,1]; % path type
% max expcected TU for hohmann transfers:
% e <-> m: 8.907 TU
% m-> j: 38.732 TU
% j-> n: 465.9 TU

% differential evolution options
opts.pop_size = 200;
opts.max_gen  = 100;
opts.F        = 0.8;
opts.CR       = 0.9;

fitness_fn = @(x) fitness(x);
[x_best, dv_best] = differential_evolution(fitness_fn, lb, ub, opts);

fprintf('Best solution: %.3f km/s total ΔV\n', dv_best);
fprintf('best state:');
disp(x_best);

% Polish with gradient descent
% options = optimset('MaxFunEvals', 50000, 'TolFun', 1e-6);
% [x_polished, dv_polished] = fminsearch(fitness_fn, x_best, options);
% 
% fprintf('true best: %.3f km/s total ΔV\n', dv_polished);
% fprintf('best state:');
% disp(x_polished);
