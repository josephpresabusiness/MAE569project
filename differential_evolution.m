function [x_best, dv_best] = differential_evolution(fitness_fn, lb, ub, opts)
    % im ngl i clauded tf out of this function
    % opts.pop_size  : e.g. 200
    % opts.max_gen   : e.g. 1000
    % opts.F         : mutation scale, e.g. 0.8
    % opts.CR        : crossover rate, e.g. 0.9

    N  = length(lb);
    NP = opts.pop_size;

    % TODO handle boolean values in state vector
    % Initialize population uniformly in [lb, ub]
    pop = lb + rand(NP, N) .* (ub - lb);

    % Evaluate initial fitness
    fitness_vals = arrayfun(@(i) fitness_fn(pop(i,:)), 1:NP);

    for gen = 1:opts.max_gen
        for i = 1:NP
            % --- Mutation: pick 3 distinct random agents (not i) ---
            candidates = setdiff(1:NP, i);
            idx = candidates(randperm(length(candidates), 3));
            a = pop(idx(1),:);  b = pop(idx(2),:);  c = pop(idx(3),:);

            mutant = a + opts.F * (b - c);
            mutant = min(max(mutant, lb), ub);   % enforce bounds

            % --- Crossover ---
            mask  = rand(1, N) < opts.CR;
            mask(randi(N)) = true;               % ensure at least 1 dimension crosses
            trial = pop(i,:);
            trial(mask) = mutant(mask);

            % --- Selection ---
            trial_fit = fitness_fn(trial);
            if trial_fit < fitness_vals(i)
                pop(i,:)        = trial;
                fitness_vals(i) = trial_fit;
            end
        end

        if mod(gen, 50) == 0
            fprintf('Gen %d | Best ΔV: %.4f km/s\n', gen, min(fitness_vals));
        end
    end

    [dv_best, best_idx] = min(fitness_vals);
    x_best = pop(best_idx,:);
end