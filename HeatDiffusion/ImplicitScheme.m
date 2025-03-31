function [T_profile, t_final] = ImplicitScheme(params)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solves 1D Heat Equation using Implicit Scheme.
    % Returns:
    %         T_profile = Final temperature profile
    %         t_final = Time at which final temperature reaches T_bondline
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Unpack parameters
    T_hot = params.T_hot;
    T_bondline = params.T_bondline;
    alpha = params.alpha;
    r = params.r;  % Stability parameter r = alpha*dt/dx^2
    Nx = params.Nx;
    dx = params.dx;
    dt = params.dt;
    max_iter = params.max_iter;

    % Initialize
    T = params.T_init * ones(Nx, 1); 
    t = 0;
    iter = 0;

    % Apply Initial Boundary Conditions
    T(1) = T_hot; % Dirichlet BC
    T(end) = T(end-1); % Neumann BC

    % Thomas Algorithm Coefficients
    a = -r * ones(Nx-1,1);          % Upper Diagonal
    d = (1+2*r) * ones(Nx,1);       % Main Diagonal
    b = -r * ones(Nx-1,1);          % Lower Diagonal
    c = T;
    g=0;                            % gradient=0

    % Modify Thomas Algorithm Coefficients for BC
    % Dirichlet at x=0
    a(1) = 0;
    d(1) = 1;
    c(1) = T_hot;

    % Neumann at x=L (zero gradient)
    b(end) = d(end-1) + 4*b(end-1);
    d(end) = a(end-1) - 3*b(end-1);
    c(end) = c(end-1) - 2*dx*b(end-1)*g;

    

    a_prime = a;
    b_prime = b;
    c_prime = c;
    d_prime = d;
    
    while T(end) < T_bondline && iter < max_iter
        iter = iter+1;

        % Forward Sweep
        % Loop over internal points
        for i = 1:Nx-1
            frac = b(i) / d_prime(i);
            d_prime(i+1) = d(i+1) - frac*a(i);
            b_prime(i) = 0;
            c_prime(i+1) = c(i+1) - frac*c_prime(i);
        end

        % Backward Sub
        % Solve for last equation
        T(Nx) = c_prime(Nx) / d_prime(Nx);
        % Loop over A matrix, but backwards
        for i = Nx-1:-1:1
            T(i) = (c_prime(i) -a(i) * T(i+1)) / d_prime(i);
        end

        % Apply BCs
        T(1) = T_hot;
        T(end)=T(end-1);

        c=T;
        c_prime = c;
        t = t+dt;
    end
    
    % Compute final time with linear interpolation
    if T(end-1) >= T_bondline
        t_final = t - dt + (T_bondline - T(end-2)) / (T(end-1) - T(end-2)) * dt;
    else
        t_final = NaN;
        error('Simulation reached max iterations before bondline temperature was achieved.');
    end
    fprintf('Time to %.2f K: %.4f seconds for Nx: %d\n', T_bondline, t_final, Nx);
    T_profile = T;
end     %end of function