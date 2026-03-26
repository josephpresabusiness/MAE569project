% TODO get some better state vector i made ts up
lb = [9497,   30,  30, 100,  500,   200,   200,  1e5,  0,  0,  0,  0];
ub = [16861, 400, 400, 800, 5000, 10000, 10000, 3e5,  1,  1,  1,  1];

% differential evolution options
opts.pop_size = 200;
opts.max_gen  = 1000;
opts.F        = 0.8;
opts.CR       = 0.9;

fitness_fn = @(x) fitness(x);
[x_best, dv_best] = differential_evolution(fitness_fn, lb, ub, opts);

fprintf('Best solution: %.3f km/s total ΔV\n', dv_best);

% Polish with gradient descent
options = optimset('MaxFunEvals', 50000, 'TolFun', 1e-6);
[x_polished, dv_polished] = fminsearch(fitness_fn, x_best, options);

fprintf('true best: %.3f km/s total ΔV\n', dv_polished);
fprintf('best state:');
disp(x_polished);