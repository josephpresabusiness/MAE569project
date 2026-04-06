function [x_best, dv_best,x_best_list,dv_best_list] = differential_evolution(fitness_fn, lb, ub, opts)
    % opts.pop_size  : e.g. 200
    % opts.max_gen   : e.g. 1000
    % opts.F         : mutation scale, e.g. 0.8
    % opts.CR        : crossover rate, e.g. 0.9

    N  = length(lb);
    NP = opts.pop_size;

    % Initialize population uniformly in [lb, ub]
    pop = lb + rand(NP, N) .* (ub - lb);

    % last four columns can only br 0,1, so rounding them will give a
    % random distribution that still works properly. 
    pop(:, end-3:end) = round(pop(:, end-3:end));

    % Evaluate initial fitness
    fitness_vals = arrayfun(@(i) fitness_fn(pop(i,:)), 1:NP);

    for gen = 1:opts.max_gen
        for i = 1:NP
            candidates = setdiff(1:NP, i);
            idx = candidates(randperm(length(candidates), 3));
            a = pop(idx(1),:);  b = pop(idx(2),:);  c = pop(idx(3),:);

            mutant = a + opts.F * (b - c);
            mutant = min(max(mutant, lb), ub); % enforce bounds

            mask  = rand(1, N) < opts.CR;
            mask(randi(N)) = true; % ensure at least 1 dimension crosses
            trial = pop(i,:);
            trial(mask) = mutant(mask);

            trial_fit = fitness_fn(trial);
            if trial_fit < fitness_vals(i)
                pop(i,:)        = trial;
                fitness_vals(i) = trial_fit;
            end
        end
        fitness_vals = fitness_vals.*29.79;
        [dv_current,idx_current] = min(fitness_vals);
        x_current = pop(idx_current,:);
        x_best_list(1,gen) = gen;
        dv_best_list(1,gen) = gen;
        dv_best_list(2,gen) = dv_current;
        x_best_list(2:13,gen) = x_current';
        fprintf('Gen %d | Best ΔV: %.4f km/s\n', gen, min(fitness_vals));
        
    end

    [dv_best, best_idx] = min(fitness_vals);
    x_best = pop(best_idx,:);
end
