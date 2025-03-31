clear;clc;close all;

%% Heat Diffusion Simulation
% This code is used to predict time for heat of re-entry to pass through
% shuttle tile and heat bondline to 600F using 1D heat equation using:
%
% a) 2-stage Runge Kutte Explicit scheme
% b) Simple Implicit numerical scheme
%
% The heat equation is given as:
% dT/dt = alpha * d^2/dx^2(T), where alpha = k/(rho*c_v)
params = struct();
%% Main Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physics Modeling
params.rho = 140;              % Density of Tile kg/m^3
params.c = 628;                % Specific Heat J/kg-K
params.k = 0.048;              % Thermal Conductivity W/m-K
params.thickness_inches = 2;   % Tile Thickness in inches
params.T_hot = 1530;           % Temperature on Re-entry Side in Kelvin
params.T_init = 300;           % Temperature at t=0 in Kelvin
params.T_bondline = 587;       % Temperature at Bondline in Kelvin
scheme = 'RK2';
% scheme = 'Implicit';

% Solver Settings
% Stability Criterion
if strcmp(scheme,'RK2')
    fprintf('Starting 2-Stage Runge Kutta Explicit Simulations\n')
    Nx_values = [2000,3000,3500];        % Number of Divisions
    params.r = 0.5;
elseif strcmp(scheme, 'Implicit')
    fprintf('Starting Simple Implicit Simulations\n')
    Nx_values = [200,500,1200];        % Number of Divisions
    params.r=2.0;
end
params.max_iter = 999999999;       % Max Iterations for debugging
T_profile_Results = cell(1,length(Nx_values));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Calculations
for i=1:length(Nx_values)
    params.Nx = Nx_values(i);
    params.L = params.thickness_inches/39.37;           % Tile Thickness in meters
    params.alpha = params.k/(params.rho*params.c);      % Thermal Diffusivity 
    params.dx = params.L / double(params.Nx-1);               % Grid spacing in meters
    params.dt = params.r*params.dx^2 / params.alpha;    % Timestep size
    params.x = linspace(0,params.L,params.Nx);          % Grid Points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(scheme,'RK2')
        %% 2-Stage Runge Kutta Explicit Scheme
        [T_final, t_results(i)] = RK2Solve(params);
    elseif strcmp(scheme,'Implicit')
        [T_final, t_results(i)] = ImplicitScheme(params);
    end

end
GridStudyResults = table(Nx_values', t_results', 'VariableNames', {'Grid_Size', 'Time_to_Bondline'});
disp(GridStudyResults);

figure;
plot(params.x, T_final, 'LineWidth', 2);
xlabel('Tile Thickness (m)');
ylabel('Temperature (K)');
titlename = sprintf('Temperature Distribution at Final Time: %.2f s for %s scheme', t_results(end), scheme);
title(titlename);
grid on;







