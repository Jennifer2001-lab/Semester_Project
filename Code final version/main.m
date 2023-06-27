close all;
clear;
Sb = 600000; Vb = 400; lpsize = [7,13]; lines = 'linedata_microgrid_original.txt';
createGridData(lines,Sb,Vb,lpsize)
createResourceParameters()
load MicroGrid.mat
load Resourceparameters.mat 
load microloadsreal.mat load_indices uncon_indices
load Season_scenarios_weeks_Jennifer.mat
load Load_scenarios_weeks_Jennifer.mat 

% addpath(genpath('C:\Users\admin\Documents\rgupta\Yalmip\'))
% addpath('C:\Program Files\Mosek\10.0\toolbox\r2017a')
% addpath('C:\gurobi1001\win64\matlab\')

%% Load the data
loads.data = load_scenarios;
loads.indices = load_indices;
uncontrollable.data = pv_scenarios;
uncontrollable.indices = uncon_indices;

el_price = elpricenextday; 
Sb = 600000;

T=7*24; % number of timesteps
S=5; % number of scenarios, max 20
subprob_number=2; % number of subproblems, max 4

% Benders decomposition settings
max_iterations = 1000;

% Initialize data storage structures
duals_all = cell(1, max_iterations);
subprob_sol_all = cell(1, max_iterations);
subprob_objective_all = cell(1, max_iterations);
master_objective_all = zeros(1, max_iterations);
alpha_all = cell(1, max_iterations);
z_up_all = zeros(1, max_iterations);
z_down_all = zeros(1, max_iterations);
eps = zeros(subprob_number, max_iterations);

for iter = 1:max_iterations
    clear gurobi; yalmip clear;
    disp(['Iteration: ', num2str(iter)]);
    
    % Solve master problem
    [master_objective, master_solution, master_alpha]  = solve_master_problem_sequence(duals_all, subprob_sol_all, subprob_objective_all, iter, subprob_number);

    new_subprob_objective = zeros(1,subprob_number);
    new_subprob_duals = cell(1,subprob_number);
       
%     for k = 1:subprob_number
%             [subproblem_objective, subprob_solution_temp,  dual_solution_sub] = solve_subproblem(loads, uncontrollable, grid, resources, costs, el_price, master_solution, T, S, k, k);
%             new_subprob_objective(k) = subproblem_objective;
%             new_subprob_duals{k} = dual_solution_sub;
%             eps_temp(k) = subprob_solution_temp.eps;
%     end

    for k = 1:subprob_number
        [subproblem_objective, subprob_solution_temp, dual_solution_sub] = solve_subproblem_sequence_matrix(loads, uncontrollable, grid, resources, costs, el_price, master_solution, T, S, k, k, k,subprob_number);
        new_subprob_objective(k) = subproblem_objective;
        new_subprob_duals{k} = dual_solution_sub;
        eps_temp(k) = subprob_solution_temp.eps;
    end

    fprintf('Ebmax: %.6f\n', master_solution.Ebmax*Sb);
    fprintf('Sbmax: %.6f\n', master_solution.Sbmax*Sb);
    fprintf('Ehmax: %.6f\n', master_solution.Ehmax*Sb);
    fprintf('dEh: %.6f\n', master_solution.dEh*Sb);
    fprintf('dEhp: %.6f\n', master_solution.dEhp*Sb);
    fprintf('dEhm: %.6f\n', master_solution.dEhm*Sb);
    
    fprintf('PFCmax: %.6f\n', master_solution.PFCmax*Sb);
    fprintf('PELmax: %.6f\n', master_solution.PELmax*Sb);


    % Store the data at the current iteration
    subprob_sol_all{iter} = master_solution;
    duals_all{iter} = new_subprob_duals;
    subprob_objective_all{iter} = new_subprob_objective;
    master_objective_all(iter) = master_objective;
    alpha_all{iter} = master_alpha;
    eps(:, iter) = eps_temp';
    
    % Check convergence
    z_up = master_objective - sum(master_alpha) + sum(new_subprob_objective); %master_obj + sub_obj
    z_down = master_objective; %master_obj + alpha

    % Store z_up and z_down values
    z_up_all(iter) = z_up;
    z_down_all(iter) = z_down;
    
    fprintf('z_up: %.3f\n', value(z_up));
    fprintf('z_down: %.3f\n', value(z_down));
   
    if abs(z_up - z_down)/abs(z_down) < 1e-5
        disp('Convergence reached');
        break;
    end
    
    if z_up<z_down
        disp('Something is wrong :( ');
        break;
    end 
       
    
    close all

    % Plot z_up and z_down

%     start=1;
%     if iter > 20
%         start = 20;
%     end

%     figure;
%     plot(start:iter, z_up_all(start:iter), 'b-', 'LineWidth', 2);
%     hold on;
%     plot(start:iter, z_down_all(start:iter), 'r--', 'LineWidth', 2);
%     xlabel('Iteration');
%     ylabel('Objective value');
%     legend('z_{up}', 'z_{down}');
%     title('Convergence of z_{up} and z_{down}');
%     grid on;


end

% Plot z_up and z_down
start = 4;
figure;
plot(start:iter-1, z_up_all(start:iter-1), 'b-', 'LineWidth', 2);
hold on;
plot(start:iter-1, z_down_all(start:iter-1), 'r--', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Objective value');
legend('z_{up}', 'z_{down}');
title('Convergence of z_{up} and z_{down}');
grid on;

%plot slack variable variation
figure;
plot(2:iter-1, eps(2:iter-1), 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Slack');
title('Slack variable variation');
grid on;
 
Ebmax_values = zeros(1, iter);
Sbmax_values = zeros(1, iter);
Ehmax_values = zeros(1, iter);

Pfc_max_values = zeros(1, iter);
Pel_max_values = zeros(1, iter);

for i = 1:iter-1
    Ebmax_values(i) = subprob_sol_all{i}.Ebmax;
    Sbmax_values(i) = subprob_sol_all{i}.Sbmax;
    Ehmax_values(i) = subprob_sol_all{i}.Ehmax;
    Pfc_max_values(i) = subprob_sol_all{i}.PFCmax;
    Pel_max_values(i) = subprob_sol_all{i}.PELmax;
end

figure
subplot(2, 3, 1)
plot(1:iter, Ebmax_values*Sb/1000, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Ebmax [kWh]');
title('Ebmax variation');
grid on;
subplot(2, 3, 2)
plot(1:iter, Sbmax_values*Sb/1000, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Sbmax [kW]');
title('Sbmax variation');
grid on;
subplot(2, 3, 4)
plot(1:iter, Ehmax_values*Sb/1000, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Ehmax [kWh]');
title('Ehmax variation');
grid on;
subplot(2, 3, 5)
plot(1:iter, Pfc_max_values*Sb/1000, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('PFCmax [kW]');
title('PFCmax variation');
grid on;
subplot(2, 3, 6)
plot(1:iter, Pel_max_values*Sb/1000, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('PELmax [kW]');
title('PELmax variation');
grid on;

fprintf('Ebmax: %.6f\n', master_solution.Ebmax*Sb);
fprintf('Sbmax: %.6f\n', master_solution.Sbmax*Sb);
fprintf('Ehmax: %.6f\n', master_solution.Ehmax*Sb);

fprintf('PFCmax: %.6f\n', master_solution.PFCmax*Sb);
fprintf('PELmax: %.6f\n', master_solution.PELmax*Sb);



%% Plot loads and seasons

plot_loads_seasons = 1;
if plot_loads_seasons == 1 
    S=15;
    time = [linspace(0,24*7,24*7)];
 
    for i=1:4
        figure
        plot(time,10e3*uncontrollable.data{i}(:,1:S)/1000) %all scenarios
        xlabel('Time (hours)')
        ylabel('Power [kW]')
        title(['PV ', num2str(i)])
        grid on
    end

    figure 
    plot(time,10e3*loads.data(:,1:S)/1000)
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    title('Load 1')
    grid on

    figure 
    plot(time,5e3*loads.data(:,1:S)/1000)
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    title('Load 2')
    grid on

    figure 
    plot(linspace(1, 24, 24), el_price(1:24))
    xlabel('Time (hours)')
    ylabel('Price [CHF/kWh]')
    title('Electricity prices')
    grid on

end 

%%
sol = subprob_solution_temp;
scale = Sb/1000;
figure(11)
hold on
%plot(sol.dispatch*scale)
plot(sol.ps*scale)
title('Dispatch Tracking')
xlabel('Time(step)')
ylabel('Power [kW]')
legend('Dispatch', 'Slack Power of Scenarios')

figure(12)
xlabel('Time(step)')
subplot(3,1,1)
ylabel('Power [kW]')
title('Battery')
hold on
plot(sol.Pb*scale)
subplot(3,1,2)
title('Fuel Cell')
ylabel('Power [kW]')
hold on
plot(sol.Pfc*scale)
subplot(3,1,3)
title('Electrolyzer')
ylabel('Power [kW]')
hold on
plot(sol.Pel*scale)

figure(13)
subplot(2,1,1)
hold on
title('Battery Energy')
ylabel('Energy [kWh]')
plot(sol.Eb*scale)
subplot(2,1,2)
hold on
title('Hydrogen Storage')
plot(sol.Eh*scale)
xlabel('Time(step)')
ylabel('Energy [kWh]')




