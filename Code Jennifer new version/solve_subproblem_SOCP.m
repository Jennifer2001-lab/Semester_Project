function [objective, solution, dual_solution_sub] = solve_subproblem_SOCP(loads, uncontrollable, grid, resources, costs, el_price, master_sol, T, S, pv_num, subprob_number)
    %% Unpack the Grid Data
    state.SOCH = 0.5; state.SOCB = 0.5; % Initial battery level, initial hydrogen pressure

    bv =grid.basevalues; % Base values
    Vb=bv(1); Sb=bv(2); pb=bv(5); Ib = Sb/(sqrt(3)*Vb);
    LP = grid.lineparameters; % Line parameters, indices, admittances, ampacities ...
    NL = grid.lines; N = grid.nodes; NC = 5; % Number of lines and nodes, controllable: battery (p and q), fuel cell (p and q), electrolyzer (q)
    idx = grid.node_indices; % nodal indices
    adm = grid.Admittance; % Admittance matrices
    YY = adm{1}; YYL=adm{2}; YYT=adm{3}; Ib=adm{4}; 
    idxbat = 15; idbat =1;
    idxh = 13; idfc = 2; idel = 3;
    Imax = repmat((LP(:,7)/Ib),T,1); vmin = 0.95; vmax = 1.05;
    paramsNR.tol=1e-8; paramsNR.n_max=1000;  % maximum number of iterations, convergence tolerance


    batteryindices = 1:T; fcindices = batteryindices; 
    batteryindicesc = batteryindices; fcindicesc = T+1:2*T; elindices = 1:T;
    fcnext = fcindices(2:end); fcprev = fcindices(1:end-1); elnext = elindices(2:end); elprev = elindices(1:end-1); batnext = batteryindices(2:end); batprev = batteryindices(1:end-1); 
    
    E_star = ones(N,1); E_0 = E_star;% Only first relevant (slack node)

    topology.downstream_lines = grid.G; topology.upstream_nodes = grid.Up; topology.R = grid.resistances(:,1); topology.X = grid.resistances(:,2); topology.B = grid.resistances(:,3);
    topology.downstream_nodes = grid.Down;
    topology.Inj = zeros(N,1); topology.Inj(idxbat) = 1; topology.Inj(idxh) = 2;
    topology.Inj(loads.indices(1)) = 3; topology.Inj(loads.indices(2)) = 4;
    topology.Inj(uncontrollable.indices(1)) = 5; topology.Inj(uncontrollable.indices(2)) = 6;
    topology.Ninj = max(topology.Inj); Ninj = topology.Ninj;


    %% Define subproblem variables
    Dslackpos = sdpvar(T,S); 
    Dslackneg = sdpvar(T,S);
    dispatch= sdpvar(T,1);
    mhprod = sdpvar(T,S); 
    mhcons = sdpvar(T,S); 
    Eh = sdpvar(T+1,S); 
    Eb = sdpvar(T+1,S); 
    Pinj = sdpvar(T*6,S); Qinj = sdpvar(6*T,S); % Active and reactive together
    Plines = sdpvar(T*NL,S); Qlines = sdpvar(T*NL, S); flines = sdpvar(T*NL,S); Vnodes = sdpvar(T*N,S);
    Pfc = sdpvar(T,S); Pel = sdpvar(T,S); Ps_pos = sdpvar(T,S); Ps_neg = sdpvar(T,S);


    Ebmax = sdpvar(1, 1);
    Sbmax = sdpvar(1, 1);
    Ehmax = sdpvar(1, 1);

    PFCmax = sdpvar(1, 1);
    PELmax = sdpvar(1, 1);


    eps = sdpvar(1, 1);

    % Unit commitment per day
%         fc_bin = sdpvar(7, 1);
%         el_bin = sdpvar(7, 1);

    % Unit commitment per week
    fc_bin = sdpvar(1, 1);
    el_bin = sdpvar(1, 1);

   
    % Hydrogen Storage Parameters
    load('MicroGrid.mat', 'grid' )
    bv= grid.basevalues;
    Vb=bv(1); Sb=bv(2); pb=bv(5); Ib = Sb/(sqrt(3)*Vb);
    Tt = 293; % Assumed tank temperature
    Ktank_prime = (10^-5) * 4124 * Tt *Ib/pb; 

    % Add constraint to use battery for daily cycling
        midnight = 1:24:T+1;
    %% Define constraints
    cons = [];
    
    cons = [cons, (Ebmax==master_sol.Ebmax) : 'Ebmax'];
    cons = [cons, (Sbmax==master_sol.Sbmax) : 'Sbmax'];
    cons = [cons, (Ehmax==master_sol.Ehmax) : 'Ehmax'];
    cons = [cons, (PFCmax==master_sol.PFCmax) : 'PFCmax'];
    cons = [cons, (PELmax==master_sol.PELmax) : 'PELmax'];

    cons = [cons, (el_bin==master_sol.el_bin(:, subprob_number)) : 'el_bin'];
    cons = [cons, (fc_bin==master_sol.fc_bin(:, subprob_number)) : 'fc_bin'];

    cons = [cons, eps >= 0];

    cons = [cons, Dslackpos >= 0, Dslackneg >= 0, Dslackneg <= 1, Dslackpos <= 1, Ps_neg >= 0, Ps_pos >= 0, Ps_pos <= resources.Pgmax, Ps_neg <= resources.Pgmax];
    % SOCP grid constraints
    for s=1:S
        cons = [cons, createconstraintsSOCP_OPF_yalmip(topology, T, Vnodes(:,s), Plines(:,s), Qlines(:,s), flines(:,s), Pinj(:,s), Qinj(:,s))];
    end
    % Bounds
    cons = [cons, Vnodes <= vmax^2, Vnodes >= vmin^2, flines >= 0, flines <= repmat(Imax.^2,1,S)];
    cons = [cons, Plines >= resources.Pgmin, Plines <= resources.Pgmax, Qlines <= resources.Pgmax, Qlines >= resources.Pgmin];

    % Enforce injections
    cons = [cons, Pinj(2*T+1:3*T,:) == loads.data(1:T,1:S)*10000/Sb];
    cons = [cons, Pinj(3*T+1:4*T,:) == loads.data(1:T,1:S)*5000/Sb];
    cons = [cons, Pinj(4*T+1:5*T,:) == -uncontrollable.data{1}(1:T,1:S)*10000/Sb];
    cons = [cons, Pinj(5*T+1:6*T,:) == -uncontrollable.data{2}(1:T,1:S)*10000/Sb];
    cons = [cons, Qinj(2*T+1:3*T,:) == zeros(T,S)];
    cons = [cons, Qinj(3*T+1:4*T,:) == zeros(T,S)];
    cons = [cons, Qinj(4*T+1:5*T,:) == zeros(T,S)];
    cons = [cons, Qinj(5*T+1:6*T,:) == zeros(T,S)];
    % Also set FC reactive power to zero
    cons = [cons, Qinj(T+1:2*T,:) == zeros(T,S)];


    [slope_E, constant_E] = Linearize_quadratic(1, 5);

    for scenario=1:S

        for k = 1:numel(slope_E)/2
            cons = [cons, -Pinj(batteryindices,scenario) <= slope_E(k)*-Qinj(batteryindicesc,scenario) + Sbmax*constant_E(k)];
            cons = [cons, -Pinj(batteryindices,scenario) >= -slope_E(k)*-Qinj(batteryindicesc,scenario) - Sbmax*constant_E(k)];
        end


          % Storage constraints
        cons = [cons, Eh(1,scenario)== state.SOCH * Ehmax];
        cons = [cons, Eb(1,scenario) == state.SOCB * Ebmax];
        cons = [cons, Eb(midnight,scenario) <= 0.55 * Ebmax];
        cons = [cons, Eb(midnight,scenario) >= 0.45 * Ebmax];
        cons = [cons, Eh(T+1,scenario) <= 0.55 * Ehmax + eps];
        cons = [cons, Eh(T+1,scenario) >= 0.5 * Ehmax - eps];
        

        M = 15; % should be bigger?
        
        % Unit commitment per day
%             for d=0:6
%                 cons = [cons, c(fcindices(1+d*24:d*24+24),scenario) <= fc_bin(1+d)*M];
%                 cons = [cons, c(fcindices(1+d*24:d*24+24),scenario) >= PFCmax*0.1 - (1 - fc_bin(1+d))*M];
%                 cons = [cons, c(elindices(1+d*24:d*24+24),scenario) >= -el_bin(1+d)*M];
%                 cons = [cons, c(elindices(1+d*24:d*24+24),scenario) <= -PELmax*0.1 + (1 - el_bin(1+d))*M]; 
%             end

        
        % Unit commitment per week
        cons = [cons, Pfc(:,scenario) <= fc_bin*M];
        cons = [cons, Pfc(:,scenario) >= PFCmax*0.1 - (1 - fc_bin)*M];
        cons = [cons, Pel(:,scenario) <= el_bin*M];
        cons = [cons, Pel(:,scenario) >= PELmax*0.1 + (1 - el_bin)*M]; 
        

        cons = [cons, Pfc(:,scenario) <= PFCmax];
        cons = [cons, Pfc(:,scenario) >= 0];
        cons = [cons, Pel(:,scenario) <= PELmax];
        cons = [cons, Pel(:,scenario) >= 0];

        % Power balance at the hydrogen node (2 injections)
        cons = [cons, Pel - Pfc == Pinj(T+1:2*T,:)];

        cons = [cons, resources.kh_fc * Pfc(:,scenario) == mhcons(:,scenario), resources.kh_el * Pel(:,scenario) == mhprod(:,scenario)];

        % Battery State of Charge
        cons = [cons, Eb(2:end,scenario) == Eb(1:end-1,scenario) - Pinj(1:T,scenario) ];% Battery energy level (assuming timesteps of 15 minutes)
        cons = [cons, Eb(:,scenario) >= Ebmax*0.2, Eb(:,scenario) <= Ebmax]; %Storage must remain above minimum level and below max

        %cons = [cons, c(qindices,:) == 0];
        % Hydrogen Storage
        cons = [cons, Eh(2:end,scenario) == Eh(1:end-1,scenario) -  Ktank_prime * (mhcons(:,scenario) - mhprod(:,scenario))]; % Hydrogen tank pressure
        cons = [cons, Eh(:,scenario) >= (1/15)*Ehmax - eps, Eh(:,scenario) <= Ehmax + eps]; % Hydrogen tank pressure limits


        % Ramping constraints
        cons = [cons, Pfc(2:T, scenario) - Pfc(1:T-1, scenario) <= 0.5 * PFCmax];
        cons = [cons, Pfc(1:T-1, scenario) -Pfc(2:T, scenario) <= 0.5 * PFCmax];
        
        cons = [cons, Pel(2:T, scenario) - Pel(1:T-1, scenario) <= 0.5 * PELmax];
        cons = [cons, Pel(1:T-1, scenario) - Pel(2:T, scenario) <= 0.5 * PELmax];
%Increased as timesteps are now hourly      
        % Slack Power and Limit on Power Factor
        cons = [cons, Plines(1:T,scenario) == Ps_pos(:,scenario) - Ps_neg(:,scenario)];
        cons = [cons, Plines(1:T,scenario) == dispatch + Dslackpos(:,scenario) - Dslackneg(:,scenario)];
        cons = [cons, Ps_pos(:,scenario) + Ps_neg(:,scenario) >= Qlines(1:T,scenario) * tan(pi/2-resources.theta_max), Ps_pos(:,scenario) + Ps_neg(:,scenario) >= - Qlines(1:T,scenario) * tan(pi/2-resources.theta_max)]; % Power Factor Constraint
    
    end

      %% Write Objective
    obj=0;

    for scenario=1:S
        
        obj = obj + 10*6* sum(max(el_price) * (Dslackpos(:,scenario) + Dslackneg(:,scenario)));

        obj = obj + 6 * sum(el_price * Ps_pos(:,scenario) - 1/3* el_price *  Ps_neg(:,scenario)); %Sb/1e6

        obj = obj + costs.cfc * (sum(Pfc(:, scenario).^2)) + costs.cel * sum(Pel(:, scenario).^2) + ...
            costs.cbat * (sum(Pinj(1:T, scenario).^2 + 100*Qinj(1:T, scenario).^2));
       
    end
    obj = obj * 1/S; % Divide by # scenarios, assuming all have the same probability
    obj = obj + 100 * sum(sum(flines)) * 1/S;
    
     obj = obj + eps*1e7;
    

    obj = obj * (10*365)/8; 


    %% Solve 

    %disp('Solving Optimization Problem (Subproblem)')
    options = sdpsettings('solver','gurobi','gurobi.Method',2,'verbose',1,'debug',1, 'gurobi.QCPDual', 1); %,'gurobi.BarQCPConvTol',1e-5);
    %options = sdpsettings('solver','mosek','verbose',0);
    sol = optimize(cons, obj, options);
    prob = sol.problem;

    dual_solution_sub.Ebmax = dual(cons('Ebmax'));
    dual_solution_sub.Sbmax = dual(cons('Sbmax'));
    dual_solution_sub.Ehmax = dual(cons('Ehmax'));
    dual_solution_sub.PFCmax = dual(cons('PFCmax'));
    dual_solution_sub.PELmax = dual(cons('PELmax'));

    dual_solution_sub.fc_bin = dual(cons('fc_bin'));
    dual_solution_sub.el_bin = dual(cons('el_bin'));


    dual_solution_sub.warning = prob;

    % Check if optimization was successful
    if prob == 0 
        %|| prob == 4
        objective = value(obj); 

        solution.Ebmax=value(Ebmax);
        solution.Sbmax=value(Sbmax);
        solution.Ehmax=value(Ehmax);
        solution.PFCmax=value(PFCmax);
        solution.PELmax=value(PELmax);

        solution.el_bin=value(el_bin);
        solution.fc_bin=value(fc_bin);

        solution.eps = value(eps);

    else
        % Infeasible or unbounded
        objective = 0;

        solution = master_sol;

        disp('Error in the optimization');

    end
    


    %% Check Loadflow
    for scenario=1:S
        ErrorsS = zeros(T,S); psbar = zeros(T,S); Vbar = zeros(T*N,S);
        E_0 = ones(N,1);
        sol.ps(:,scenario) = value(Plines(1:T,scenario)); 
        sol.qs(:,scenario) = value(Qlines(1:T,scenario)); 
        sol.Pb(:,scenario) = -value(Pinj(1:T,scenario));
        sol.Qb(:,scenario) = -value(Qinj(1:T,scenario));

        sol.Pfc(:,scenario) = value(Pfc(:,scenario)); 
        sol.Qfc(:,scenario) = -value(Qinj(T+1:2*T,scenario)); 
        sol.Pel(:,scenario) =  value(Pel(:,scenario)); 

        sol.V = sqrt(value(Vnodes));
        sol.I = sqrt(value(flines));
        sol.mh = value(mhcons - mhprod); 
        sol.Eb = value(Eb);
        sol.Eh = value(Eh);
    end


    for t=1:T
        Sinj = zeros(N,S);
        for scenario=1:S
            Sinj(1,scenario) = sol.ps(t,scenario) + 1i * sol.qs(t,scenario);
            Sinj(loads.indices(1),scenario) = -10e3*loads.data(t,scenario)*1/Sb;
            Sinj(loads.indices(2),scenario) = -5e3*loads.data(t,scenario)*1/Sb;
            Sinj(uncontrollable.indices(1),scenario) = 10e3*uncontrollable.data{pv_num}(t,scenario)*1/Sb;
            Sinj(uncontrollable.indices(2),scenario) = 10e3*uncontrollable.data{pv_num}(t,scenario)*1/Sb;


            Sinj(idxbat,scenario) = sol.Pb(t,scenario) + 1i * sol.Qb(t,scenario);
            Sinj(idxh,scenario) = sol.Pfc(t,scenario)-sol.Pel(t,scenario) + 1i * sol.Qfc(t,scenario);
            Vt = sol.V(t:T:(N-1)*T+t,scenario);
            [J,E, Snew,n_iter] = NR_polar_sol(Sinj(:,scenario),Vt,YY,E_0,idx,paramsNR);
            psbar(t,scenario) = real(Snew(1)); qsbar(t,scenario) = imag(Snew(1));
            ErrorsS(t,scenario) = abs(psbar(t,scenario)-sol.ps(t,scenario));

            Vbar(t:T:(N-1)*T+t,scenario) = E;
        end
    end
   disp('LF Errors:  ');disp(max(max(ErrorsS)));

%% Displaying errors
disp('LF:  '); disp(max(ErrorsS));

sol.ps = psbar; sol.qs = qsbar; 
sol.V = abs(Vbar); 
sol.pb = sol.Pb;
sol.D = value(dispatch);

iterationcost = value(obj);

%% Plotting

plotting = 1;
time = [linspace(0,24*7,24*7)];

if plotting ==1

    figure(1)
    for i =1:14
        subplot(4,4,i)
        Vs = reshape(sol.V(:,scenario),N,[]);
        plot(time,Vs(i,:))
        title('Nodal Voltages')
    end

    
    
    figure(2)
    hold on
    Ehs = [sol.Eh(1:end-1,scenario)];
    plot(time,Ehs*pb);
    title('Hydrogen Storage')
    xlabel('Time (hours)')
    ylabel('Pressure [bar]')
    
    figure(7)
    ActiveInj = cell(3,1); ActiveInj{1} = 'Battery'; ActiveInj{2} = 'Fuel Cell'; ActiveInj{3} =  'Electrolyzer';
    subplot(3,3,1)
    plot(time,sol.Pb*Sb/1000)
    title(ActiveInj(1))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    subplot(3,3,2)
    hold on
    plot(time,sol.Pfc*Sb/1000)
    title(ActiveInj(2))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    subplot(3,3,3)
    hold on
    plot(time,sol.Pel*Sb/1000)
    title(ActiveInj(3))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
%% plot dispatch plan
figure(3)
hold on
for i = 1:min(S,9)
    subplot(3, 3, i)
    hold on
    psum = 0;
    pss = sol.ps(:,i);
    psum=psum + pss*Sb/1000;
    plot(time, pss*Sb/1000)
    
    plot(time, value(dispatch)*Sb/1000)
    hold off
    title(['Slack Power - Scenario ', num2str(i)])
    xlabel('Time (hours)')
    ylabel('Power [kW]')
end


figure(4)
for s=1:S
    subplot(3,3,s)
    hold on
    plot(psbar(:,s))
    plot(sol.ps(:,s))
    title('Slack Power')
    legend('SOCP', 'LF')
    xlabel('Time (hours)')
    ylabel('Power [kW]')
end

figure(5)
hold on
plot(sol.D)
plot(sol.ps)
title('Dispatch Tracking')
legend('Dispatch', 'Slack')
xlabel('Time (hours)')
ylabel('Power [kW]')

figure(6)
hold on
plot(sol.D)
plot(psbar)
title('Dispatch Tracking (LF Solution)')
legend('Dispatch', 'Slack (LF)')
xlabel('Time (hours)')
ylabel('Power [kW]')


fprintf('Ebmax: %.6f\n', value(Ebmax)*Sb);
fprintf('Sbmax: %.6f\n', value(Sbmax)*Sb);
fprintf('Ehmax: %.6f\n', value(Ehmax)*Sb);

fprintf('PFCmax: %.6f\n', value(PFCmax)*Sb);
fprintf('PELmax: %.6f\n', value(PELmax)*Sb);

end