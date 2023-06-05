close all;
clear;

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

T=24*7; % number of timesteps
S=4; % number of scenarios, max 8
subprob_number=4; %number of subproblems

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
    disp(['Iteration: ', num2str(iter)]);
    
    % Solve master problem
    [master_objective, master_solution, master_alpha]  = solve_master_problem(duals_all, subprob_sol_all, subprob_objective_all, iter, subprob_number);

    new_subprob_objective = zeros(1,subprob_number);
    new_subprob_duals = cell(1,subprob_number);
    
    
%     for k = 1:subprob_number
%             [subproblem_objective, ~,  dual_solution_sub] = solve_subproblem(loads, uncontrollable, grid, resources, costs, el_price, master_solution, T, S, k, k);
%             new_subprob_objective(k) = subproblem_objective;
%             new_subprob_duals{k} = dual_solution_sub;
%     end

    parfor k = 1:subprob_number
        [subproblem_objective, subprob_solution_temp, dual_solution_sub] = solve_subproblem_(loads, uncontrollable, grid, resources, costs, el_price, master_solution, T, S, k, k);
        new_subprob_objective(k) = subproblem_objective;
        new_subprob_duals{k} = dual_solution_sub;
        eps_temp(k) = subprob_solution_temp.eps;
    end

   
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
start = 1;
figure;
plot(start:iter, z_up_all(start:iter), 'b-', 'LineWidth', 2);
hold on;
plot(start:iter, z_down_all(start:iter), 'r--', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Objective value');
legend('z_{up}', 'z_{down}');
title('Convergence of z_{up} and z_{down}');
grid on;

figure;
plot(start:iter, eps(start:iter), 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Slack');
title('Slack variable variation');
grid on;

%% Plot loads and seasons

plot_loads_seasons = 0;
if plot_loads_seasons == 1 
    S = 8;
    time = [linspace(0,24*7,24*7)];
    
    figure
    for i=1:4
        subplot(2,2,i)
        plot(time,10e3*uncontrollable.data{i}(:,1:S)/1000) %all scenarios
        xlabel('Time (hours)')
        ylabel('Power [kW]')
        title(['PV ', num2str(i)])
    end

    figure 
    plot(time,10e3*loads.data(:,1:S)/1000)
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    title('Load')

    figure 
    plot(linspace(1, 24, 24), el_price(1:24))
    xlabel('Time (hours)')
    ylabel('Price [CHF/kWh]')
    title('Electricity prices')

end 


function [objective, sol, alpha_] = solve_master_problem(duals_all, subprob_sol_all, subprob_objective_all, current_iter, subprob_number)
    Sb = 600000; 
    Vb = 400; 

    T = 24*4;

    %% Declare optimization variables
    alpha = sdpvar(1, subprob_number);
    
    Ebmax = sdpvar(1, 1);
    Sbmax = sdpvar(1, 1);
    Ehmax = sdpvar(1, 1);
    
    PFCmax = sdpvar(1, 1);
    PELmax = sdpvar(1, 1);
    
    hyd_install = binvar(1, 1);
    bat_install = binvar(1, 1);

    fc_bin = binvar(1, subprob_number);
    el_bin = binvar(1, subprob_number);

    %% Objective function
    obj = 0;
    obj = obj + (bat_install + hyd_install)*20000;
    obj = obj + (Ebmax + Sbmax)*400*Sb/1000 + Ehmax*300*Sb/1000 + (PFCmax + PELmax)*1000*Sb/1000;
    obj = obj + sum(alpha);

    obj = obj + 1000*sum(el_bin + fc_bin);

    %% Constraints
    cons = [];

    %initialization of alpha
    cons = [cons, alpha >= -1e9];
    
    %constraints on size 
    cons = [cons, Ebmax>=0, Sbmax>=0, Ehmax>=0, PFCmax>=0, PELmax>=0];
    cons = [cons, Sbmax<=2*Ebmax];
    
    M = 1;
    cons = [cons, Ebmax<=bat_install*M, Sbmax<=bat_install*M];
    cons = [cons, Ehmax<=hyd_install*10*M, PFCmax<=hyd_install*M, PELmax<=hyd_install*M];
    
%     cons = [cons, Ehmax <= PFCmax*20 ; Ehmax >= PFCmax*1];
%     cons = [cons, Ehmax <= PELmax*20 ; Ehmax >= PELmax*1];

    cons = [cons, PELmax <= 0.07, PFCmax <= 0.07, Ebmax <= 4*0.07];

    if current_iter > 1
        % Add Benders cuts for all saved iterations
        for i = 1:current_iter-1
            for j = 1:subprob_number
                if duals_all{i}{j}.warning == 0 
                    val_Ebmax = duals_all{i}{j}.Ebmax * (Ebmax - subprob_sol_all{i}.Ebmax);
                    val_Sbmax = duals_all{i}{j}.Sbmax * (Sbmax - subprob_sol_all{i}.Sbmax);
                    val_Ehmax = duals_all{i}{j}.Ehmax * (Ehmax - subprob_sol_all{i}.Ehmax);
                    val_PFCmax = duals_all{i}{j}.PFCmax * (PFCmax - subprob_sol_all{i}.PFCmax);
                    val_PELmax = duals_all{i}{j}.PELmax * (PELmax - subprob_sol_all{i}.PELmax);
                    
                    val_el_bin = duals_all{i}{j}.el_bin * (el_bin(j) - subprob_sol_all{i}.el_bin(j));
                    val_fc_bin = duals_all{i}{j}.fc_bin * (fc_bin(j) - subprob_sol_all{i}.fc_bin(j));

                    cons = [cons, alpha(j) >= subprob_objective_all{i}(j) - val_Ebmax - val_Sbmax - val_Ehmax - val_PFCmax - val_PELmax - val_fc_bin - val_el_bin];
          
%                     cons = [cons, alpha(j) >= subprob_objective_all{i}(j) - val_Ebmax - val_Sbmax - val_Ehmax - val_PFCmax - val_PELmax];

                else %unfeasiblity cuts
                    disp('Unfeasiblity cuts')
                    val_Ebmax = duals_all{i}{j}.Ebmax * (Ebmax - subprob_sol_all{i}.Ebmax);
                    val_Sbmax = duals_all{i}{j}.Sbmax * (Sbmax - subprob_sol_all{i}.Sbmax);
                    val_Ehmax = duals_all{i}{j}.Ehmax * (Ehmax - subprob_sol_all{i}.Ehmax);
                    val_PFCmax = duals_all{i}{j}.PFCmax * (PFCmax - subprob_sol_all{i}.PFCmax);
                    val_PELmax = duals_all{i}{j}.PELmax * (PELmax - subprob_sol_all{i}.PELmax);
                    
                    cons = [cons, - val_Ebmax - val_Sbmax - val_Ehmax - val_PFCmax - val_PELmax <= 0];
                end
            end
        end 
    end
 
    %% Solve
    % Solver settings
    options = sdpsettings('solver', 'gurobi', 'gurobi.Method', 2, 'verbose', 1, 'debug', 1, 'gurobi.QCPDual', 1);
%     options = sdpsettings('solver','mosek','verbose',0);
    % Solve the optimization problem
    solution = optimize(cons, obj, options);

    % Check for errors
    if solution.problem ~= 0
        disp('!!!!!!!!!!! Error Code !!!!!!!!!!!!');
        disp(solution.problem);
        disp('Breakpoint');
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    end
    
    
    objective = value(obj);
    
    alpha_ = value(alpha);
   
    
    sol.Ebmax = value(Ebmax);
    sol.Sbmax = value(Sbmax);
    sol.Ehmax = value(Ehmax);
    sol.PFCmax = value(PFCmax);
    sol.PELmax = value(PELmax);

    sol.el_bin = value(el_bin);
    sol.fc_bin = value(fc_bin);

    sol.battery_install = value(bat_install);
    sol.hydrogen_install = value(hyd_install);


    
end

