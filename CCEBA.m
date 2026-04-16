% Main programs starts here
function [best, Convergence_curve] = CCEBA(N, MaxFEs, lb, ub, dim, fobj)
    A = rand(1, N) + ones(1, N);  % Loudness (constant or decreasing)
    r = rand(1, N);  % Pulse rate (constant or decreasing)
    Qmin = 0;  % Frequency minimum
    Qmax = 2;  % Frequency maximum
    d = dim;  % Number of dimensions 
    Lb = lb * ones(1, d);  % Lower limit/bounds/ a vector
    Ub = ub * ones(1, d);  % Upper limit/bounds/ a vector
    Q = zeros(N, 1);  % Frequency
    v = zeros(N, d);  % Velocities
    X = zeros(N, d);  % Positions of the bats
    Fitness = zeros(N, 1);

    % Initialize the population/solutions
    for i = 1:N
        X(i,:) = Lb + (Ub - Lb) .* rand(1, d);
        Fitness(i) = fobj(X(i,:));
    end

    % Find the initial best solution
    [fmin, I] = min(Fitness);
    best = X(I, :);
    Convergence_curve = [];
    FEs = 0;
    t = 1;

    % Main loop
    while FEs < MaxFEs
                % Apply crossover mechanism
        [sorted_AllFitness, ~] = sort(Fitness);
        Rolette_index = RouletteWheelSelection(1./sorted_AllFitness);
        if Rolette_index == -1  
            Rolette_index = 1;
        end

        history_X = X;
        for i = 1:N
            Q(i) = Qmin + (Qmin - Qmax) * rand;
            v(i, :) = v(i, :) + (X(i, :) - best) * Q(i);
            S(i, :) = X(i, :) + v(i, :);

            % Apply simple bounds/limits
            S(i, :) = simplebounds(S(i, :), Lb, Ub);

            % Pulse rate
            if rand > r(i)
                % The factor 0.001 limits the step sizes of random walks 
                S(i, :) = best + A(i) * randn(1, d);
            end

            % Apply simple bounds/limits
            S(i, :) = simplebounds(S(i, :), Lb, Ub);

            if FEs < MaxFEs
                FEs = FEs + 1;
                Fnew = fobj(S(i, :));
                % Update if the solution improves, or not too loud
                if Fnew <= Fitness(i)
                    X(i, :) = S(i, :);
                    Fitness(i) = Fnew;
                end
                if Fnew < fmin && rand < A(i)
                    best = S(i, :);
                    fmin = Fnew;
                end
            else
                break;
            end
            
        end

        % Crossover
        for i = 1:N
            m = zeros(1, d);
            u = randperm(d);
            m(1, u(1:floor(rand * d + 1))) = 1;
            X2(i, :) = history_X(i, :) + randn * m .* (history_X(Rolette_index, :) - history_X(i, :));
        end

        for i = 1:N
            % Return back the search agents that go beyond the boundaries of the search space
            Flag4ub = X(i,:) > Ub;
            Flag4lb = X(i,:) < Lb;
            X(i,:) = (X(i,:) .* (~(Flag4ub + Flag4lb))) + Ub .* Flag4ub + Lb .* Flag4lb;

            Flag4ub = X2(i,:) > Ub;
            Flag4lb = X2(i,:) < Lb;
            X2(i,:) = (X2(i,:) .* (~(Flag4ub + Flag4lb))) + Ub .* Flag4ub + Lb .* Flag4lb;

            % Calculate objective function for each search agent
            Fitness(i) = fobj(X(i,:));
            FEs = FEs + 1;

            fitness1 = fobj(X2(i, :));
            FEs = FEs + 1;

            if fitness1 < Fitness(i)
                Fitness(i) = fitness1;
                X(i, :) = X2(i, :);
            end

            if Fitness(i) < fmin
                fmin = Fitness(i);
                best = X(i, :);
            end
        end
         [X, Fitness, fg] = CC(X, Fitness, d, lb, ub, fobj);
        FEs = FEs + fg;
        Convergence_curve(t) = fmin;
        t = t + 1;
    end
end

% Application of simple limits/bounds
function s = simplebounds(s, Lb, Ub)
    % Apply the lower bound vector
    ns_tmp = s;
    I = ns_tmp < Lb;
    ns_tmp(I) = Lb(I);

    % Apply the upper bound vector 
    J = ns_tmp > Ub;
    ns_tmp(J) = Ub(J);
    % Update this new move 
    s = ns_tmp;
end

% Helper function: Roulette Wheel Selection
function choice = RouletteWheelSelection(weights)
    accumulation = cumsum(weights);
    p = rand() * accumulation(end);
    chosen_index = -1;
    for index = 1:length(accumulation)
        if accumulation(index) > p
            chosen_index = index;
            break;
        end
    end
    choice = chosen_index;
end

% Crossover mechanism
function [X, fitness, fg] = CC(X, fitness, dim, lb, ub, fobj)
    fg = 0;
    Mhc = zeros(size(X, 1), dim);
    Bhc = randperm(size(X, 1));
    for i = 1:(size(X, 1) / 2)
        no1 = Bhc(2 * i - 1);
        no2 = Bhc(2 * i);
        for j = 1:dim
            r1 = unifrnd(0, 1);
            r2 = unifrnd(0, 1);
            c1 = (rand(1) * 2) - 1;
            c2 = (rand(1) * 2) - 1;
            Mhc(no1, j) = r1 * X(no1, j) + (1 - r1) * X(no2, j) + c1 * (X(no1, j) - X(no2, j));
            Mhc(no2, j) = r2 * X(no2, j) + (1 - r2) * X(no1, j) + c2 * (X(no2, j) - X(no1, j));
        end
    end
    for i = 1:size(X, 1)
        FU = Mhc(i, :) > ub;
        FL = Mhc(i, :) < lb;
        Mhc(i, :) = (Mhc(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
        fitness_mhc(i) = fobj(Mhc(i, :));
        fg = fg + 1;
        if fitness(i) < fitness_mhc(i)
            X(i, :) = X(i, :);
        else
            X(i, :) = Mhc(i, :);
            fitness(i) = fitness_mhc(i);
        end
    end
    Bvc = randperm(dim);
    Mvc = X;
    for i = 1:size(X, 1)
        Boundary_no = size(ub, 2);
        if Boundary_no == 1
            Mvc(i, :) = (Mvc(i, :) - lb) / (ub - lb);
        end
        if Boundary_no > 1
            for j = 1:dim
                ub_j = ub(j);
                lb_j = lb(j);
                Mvc(i, j) = (Mvc(i, j) - lb_j) / (ub_j - lb_j);
            end
        end
    end
    p2 = 0.6;
    for i = 1:(dim / 2)
        p = unifrnd(0, 1);
        if p < p2
            no1 = Bvc(2 * i - 1);
            no2 = Bvc(2 * i);
            for j = 1:size(X, 1)
                r = unifrnd(0, 1);
                Mvc(j, no1) = r * Mvc(j, no1) + (1 - r) * Mvc(j, no2);
            end
        end
    end
    for i = 1:size(X, 1)
        Boundary_no = size(ub, 2);
        if Boundary_no == 1
            Mvc(i, :) = Mvc(i, :) * (ub - lb) + lb;
        end
        if Boundary_no > 1
            for j = 1:dim
                ub_j = ub(j);
                lb_j = lb(j);
                Mvc(i, j) = (ub_j - lb_j) * Mvc(i, j) + lb_j;
            end
        end
        FU = Mvc(i, :) > ub;
        FL = Mvc(i, :) < lb;
        Mvc(i, :) = (Mvc(i, :) .* (~(FU + FL))) + ub .* FU + lb .* FL;
        fitness_mvc(i) = fobj(Mvc(i, :));
        fg = fg + 1;
        if fitness(i) < fitness_mvc(i)
            X(i, :) = X(i, :);
        else
            X(i, :) = Mvc(i, :);
            fitness(i) = fitness_mvc(i);
        end
    end
end