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

    % unit commitment per day
%     fc_bin = binvar(7, subprob_number);
%     el_bin = binvar(7, subprob_number);

    % unit commitment per week
    fc_bin = binvar(1, subprob_number);
    el_bin = binvar(1, subprob_number);

    %% Objective function
    obj = 0;
    obj = obj + (bat_install + hyd_install)*20000;
    obj = obj + (Ebmax + Sbmax)*400*Sb/1000 + Ehmax*100*Sb/1000 + (PFCmax + PELmax)*1000*Sb/1000;
    obj = obj + 100*sum(sum(el_bin + fc_bin));
    obj = obj + sum(alpha);

    
    %% Constraints
    cons = [];

    %initialization of alpha
    cons = [cons, alpha >= -1e13];
    
    %constraints on size 
    cons = [cons, Ebmax>=0, Sbmax>=0, Ehmax>=0, PFCmax>=0, PELmax>=0];
    cons = [cons, Sbmax<=2*Ebmax];
    
    M = 2;
    cons = [cons, Ebmax<=bat_install*M, Sbmax<=bat_install*M];
    cons = [cons, Ehmax<=hyd_install*3*M, PFCmax<=hyd_install*M, PELmax<=hyd_install*M];
    


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
                    
                    val_el_bin = (duals_all{i}{j}.el_bin)' * (el_bin(:, j) - subprob_sol_all{i}.el_bin(:, j));
                    val_fc_bin = (duals_all{i}{j}.fc_bin)' * (fc_bin(:, j) - subprob_sol_all{i}.fc_bin(:, j));

                    cons = [cons, alpha(j) >= subprob_objective_all{i}(j) - val_Ebmax - val_Sbmax - val_Ehmax - val_PFCmax - val_PELmax - val_fc_bin - val_el_bin];

                else %unfeasiblity cuts, not working ...
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
