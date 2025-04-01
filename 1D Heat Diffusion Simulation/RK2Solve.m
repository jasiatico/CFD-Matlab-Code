function [T_profile, t_final] = RK2Solve(params)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solves 1D Heat Equation using RK2 Explicit Scheme.
    % Returns:
    %         T_profile = Final temperature profile
    %         t_final = Time at which final temperature reaches T_bondline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unpack variables
    T_hot = params.T_hot;
    T_bondline = params.T_bondline;
    alpha = params.alpha;
    Nx = params.Nx;
    dx = params.dx;
    dt = params.dt;
    max_iter = params.max_iter;

    % Initialize Temperature
    T = params.T_init * ones(1, Nx); 
    Res = zeros(1, Nx); 
    t = 0;
    iter = 0;

    % Apply Initial Boundary Conditions
    T(1) = T_hot; % Dirichlet BC
    T(end) = T(end-1); % Neumann BC

    % Time Marching Loop
    while T(end) < T_bondline && iter < max_iter
        iter = iter + 1;

        % Compute Residuals (Vectorized)
        Res(2:Nx-1) = alpha * (T(3:Nx) - 2 * T(2:Nx-1) + T(1:Nx-2)) / dx^2;

        % Predictor Step (Vectorized)
        T_P = T; % Store previous state
        T_P(2:Nx-1) = T(2:Nx-1) + (dt/2) * Res(2:Nx-1);

        % Apply Boundary Conditions to T_P
        T_P(end) = T_P(end-1); % Neumann BC

        % Compute Residuals for Corrector Step using T_P
        Res(2:Nx-1) = alpha * (T_P(3:Nx) - 2 * T_P(2:Nx-1) + T_P(1:Nx-2)) / dx^2;

        % Corrector Step (Vectorized) using T_P
        T(2:Nx-1) = T(2:Nx-1) + dt * Res(2:Nx-1);

        % Apply Boundary Conditions
        T(end) = T(end-1); % Neumann BC

        % Update Time
        t = t + dt;
    end

    % Compute Final Time to Reach Bondline Temperature
    if T(end) >= T_bondline
        t_final = t - dt + (T_bondline - T_P(end)) / (T(end) - T_P(end)) * dt;
    else
        t_final = NaN;
        error('Simulation reached max iterations before bondline temperature was achieved.');
    end

    fprintf('Time to %.2f K: %.4f seconds for Nx: %d\n', T_bondline, t_final, Nx);
    T_profile = T;

end % End of function
