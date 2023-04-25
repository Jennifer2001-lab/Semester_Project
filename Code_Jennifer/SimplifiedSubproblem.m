%%% Compute a weekahead dispatch based on given FC/EL schedule and SC %%%
%%% models, only first day kept and at 15' resolution, rest at hourly %%%
function [sol, iterationcost] = SimplifiedSubproblem(loads,uncontrollable,scenario, grid, resources, costs, state, T, dispatch)
%% Unpack the Grid Data
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

converged = 0; tol = 1e-6; niter = 0; 
psinds = linspace(1,2*T-1,T); qsinds = linspace(2,2*T,T); 
batteryindices = linspace(idbat,(T-1)*NC+idbat,T); fcindices = linspace(idfc,(T-1)*NC+idfc,T); elindices = linspace(idel,(T-1)*NC+idel,T); 
batteryindicesc = linspace(idbat+3,(T-1)*NC+idbat+3,T); fcindicesc = linspace(idfc+3,(T-1)*NC+idfc+3,T);
fcnext = fcindices(2:end); fcprev = fcindices(1:end-1); elnext = elindices(2:end); elprev = elindices(1:end-1); batnext = batteryindices(2:end); batprev = batteryindices(1:end-1); 
qindices = sort([batteryindicesc, fcindicesc]);

SC.Al = cell(1,1); SC.Avm = cell(1,1); SC.Aim = cell(1,1);

E_star = ones(N,1); E_0 = E_star;% Only first relevant (slack node)
cbar = zeros(T*NC,1); Lbar = zeros(2*T,1); sbar = Lbar;
qsbar = zeros(T,1); Vbar = zeros(N*T,1); Ibar = zeros(NL*T,1); 
Sbars = zeros(T*N,1); psbar = zeros(T,1); 
for t=1:T    
    Sinj = zeros(N,1);
    Sinj(loads.indices(1)) = -loads.data{1}(t,scenario)*1/Sb;
    Sinj(loads.indices(2)) = -loads.data{2}(t,scenario)*1/Sb;
    Sinj(uncontrollable.indices(1)) = uncontrollable.data{1}(t,scenario)*1/Sb;
    Sinj(uncontrollable.indices(2)) = uncontrollable.data{2}(t,scenario)*1/Sb;
    % Compute LoadFlow solution of initial state
    [J,E,Slf,n_iter] = NR_polar_sol(Sinj,E_star,YY,E_0,idx,paramsNR);
    Sbars((t-1)*N+1:t*N) = Slf;
    Lbar((t-1)*2+1:2*t) = [sum(real(Slf));sum(imag(Slf))];
    pt = real(Slf); qt = imag(Slf); 
    sbar(2*(t-1)+1:2*(t-1) + 2) = [pt(1);qt(1)];
    psbar(t) = pt(1); qsbar(t) = qt(1); 
    Vbar((t-1)*N+1:t*N) = E; 
    Itemp = zeros(NL,1);
    for i=1:NL
        idx1 = LP(i,1); idx2 = LP(i,2);
        Itemp(i) = YYL(idx1,idx2)*(E(idx1)-E(idx2)) + YYT(idx1,idx2)*E(idx1);
       % Itemp(i+NL) = YYL(idx2,idx1)*(E(idx2)-E(idx1)) +
       % YYT(idx2,idx1)*E(idx2); Ideally also check second side
    end
    Ibar((t-1)*NL+1:t*NL) = abs(Itemp);
end
nmax = 10; % Limit the number of iterations
alphap = zeros(1,N); alphaq = zeros(1,N); % Voltage dependency of injections (not considered here)
Res_nodes_no_slack = [2:N]; slack=1; nph=1; % Parameters for SC computation, grid connection is considered as slack node
while converged == 0 && niter <= nmax
    niter = niter+1; 
    %% Compute Sensitivity Coefficients
    for t=1:T
        E = Vbar((t-1)*N+1:t*N); St = Sbars((t-1)*N+1:t*N);
        
        % Voltage Sensitivity
        [K_p,K_com_p,K_q,K_comp_q]=Coeffs_Voltage_Alpha(YY,St,E,Res_nodes_no_slack,slack,nph,alphap,alphaq,E_0);
        % Current Sensitivity
        [Icoeff, Icoeff_complex, Icoeffq,Icoeffq_complex, Currents]=Coeffs_Currents(YYL,YYT,E,K_com_p,K_comp_q,nph);
        % Write in easier format for constraints
        Aim = zeros(NL,(N-1)); % Sensitivities for one! sides of branches + wrt p and q
        for l=1:N-1
            Ic = cell2mat(Icoeff{l});
            Icq = cell2mat(Icoeffq{l});
            for i=1:NL
                idx1 = LP(i,1); idx2 = LP(i,2);
                Aim(i,l) = Ic(idx1,idx2);
                %Aim(i+NL,l) = Ic(idx2,idx1);
                Aim(i,l+N-1) = Icq(idx1,idx2);
                %Aim(i+NL,l+N-1) = Icq(idx2,idx1);
            end
        end
        % Sensitivity of the losses
        [C_rp , C_rq, C_xp, C_xq] =Coeffs_Losses(YY, E, K_com_p, K_comp_q, Res_nodes_no_slack);
        % Save SC for every timestep and every scenario
        SC.Al{t} = [C_rp(idxbat-1),C_rp(idxh-1), C_rp(idxh-1), C_rq(idxbat-1),C_rq(idxh-1); ...
            C_xp(idxbat-1), C_xp(idxh-1), C_xp(idxh-1), C_xq(idxbat-1), C_xq(idxh-1)];
        SC.Avm{t} = [K_p(:,idxbat-1),  K_p(:,idxh-1), K_p(:,idxh-1), K_q(:,idxbat-1), K_q(:,idxh-1)];
        SC.Aim{t} = [Aim(:,idxbat-1),  Aim(:,idxh-1), Aim(:,idxh-1), Aim(:,idxbat-1 + NL),  Aim(:,idxh-1 + NL)];
    end 
    Vbar = abs(Vbar);
    
    %% Define subproblem variables
    Dslackpos = sdpvar(T,1); Dslackneg = sdpvar(T,1);
    mhprod = sdpvar(T,1); mhcons = sdpvar(T,1);
    ph = sdpvar(T+1,1); 
    Eb = sdpvar(T+1,1); 
    c = sdpvar(NC*T,1); % Active and reactive together
    s = sdpvar(2*T,1); dL = sdpvar(2*T,1);
    Rbat = sdpvar(T-1,1); Rfc = sdpvar(T-1,1); Rel = sdpvar(T-1,1); % Take into account ramping (absolute values)
    dphslack = sdpvar(1,1);
    %% Define constraints
    cons = [];
    cons = [cons, Dslackpos >= 0, Dslackneg >= 0, Dslackneg <= 1, Dslackpos <= 1];
    % Impose grid constraints using SC
    % Voltage limits
    Avsparse = sparse(blkdiag(SC.Avm{:}));
    cons = [cons, Avsparse * (c - cbar) >= vmin - Vbar]; % c-cbar represents [delta_p, delta_q]
    cons = [cons, Avsparse * (c - cbar) <= vmax - Vbar];
    % Current Limits
    Aisparse = sparse(blkdiag(SC.Aim{:}));
      cons = [cons, Aisparse * (c - cbar) >= -Imax - Ibar];
      cons = [cons, Aisparse * (c - cbar) <= Imax - Ibar];
     % Compute losses based on injections
    Alsparse = sparse(blkdiag(SC.Al{:}));
    cons = [cons, dL == Alsparse * (c - cbar)];
     % Compute power at the slack node
     Acol = [1,1,1,0,0;0,0,0,1,1]; Acolcell = repmat({Acol},1,T);
     Acolsparse = sparse(blkdiag(Acolcell{:})); 
     cons = [cons, (s-sbar) == dL - Acolsparse * (c-cbar), dL<=1, dL>= -1, s >= -1, s <= 1];
     % Slack node
     cons = [cons, s(psinds) <= resources.Pgmax, s(psinds) >= resources.Pgmin];

     % Storage constraints
     cons = [cons, ph(1)== state.SOCH * resources.phmax]; % Storage constraints 
     cons = [cons, Eb(1) == state.SOCB * resources.Ebmax];
     % Resource constraints
     cons = [cons, c(fcindices) <= resources.PFCmax, c(fcindices) >= resources.PFCmin];
     cons = [cons, c(elindices) >= -resources.PELmax, c(elindices) <= -resources.PELmin];
     cons = [cons, c(fcindicesc) <= resources.QFCmax, c(fcindicesc) >= resources.QFCmin]; 
     cons = [cons, c(batteryindices) <= resources.Pbmax, c(batteryindicesc) <= resources.Qbmax  ,...
          resources.Pbmin <= c(batteryindices), resources.Qbmin <= c(batteryindicesc)]; 
     cons = [cons, resources.kh_fc * c(fcindices) == mhcons, -resources.kh_el * c(elindices) == mhprod];
     % Batery State of Charge
     cons = [cons, Eb(2:end) == Eb(1:end-1) - c(batteryindices(1:T)) *1/4];% Battery energy level (assuming timesteps of 15 minutes)
     cons = [cons, Eb >= resources.Ebmin, Eb <= resources.Ebmax]; %Storage must remain above minimum level and below max
     % Hydrogen Storage
     cons = [cons, ph(2:end) == ph(1:end-1) - resources.Ktank * (mhcons - mhprod)*1/4]; % Hydrogen tank pressure
     cons = [cons, ph >= resources.phmin - dphslack, ph <= resources.phmax + dphslack]; % Hydrogen tank pressure limits

     cons = [cons, Rbat >= 0, Rbat >= c(batnext) - c(batprev),...
          Rbat >= c(batprev) - c(batnext)];
     cons = [cons, Rfc >= 0, Rfc >= c(fcnext) - c(fcprev),...
          Rfc >= c(fcprev) - c(fcnext)];
     cons = [cons, Rel>=0, Rel >= c(elnext) - c(elprev),...
          Rel >= c(elprev) - c(elnext)];
     
      cons = [cons, dphslack >=0, dphslack <= 1];
      % Slack Power and Limit on Power Factor
      cons = [cons, s(psinds) == dispatch + Dslackpos - Dslackneg];
      cons = [cons, abs(dispatch) >= s(qsinds) * tan(pi/2-resources.theta_max), abs(dispatch) >= - s(qsinds) * tan(pi/2-resources.theta_max)]; % Power Factor Constraint
    
    %% Write Objective
    obj = 0;
    obj = obj + dphslack *  costs.cslack;
    obj = obj +  1e5 * (sum(Dslackpos.^2) + sum(Dslackneg.^2));  % cost for deviating from dispatch %%      
    obj = obj + costs.cfc * sum(Rfc) + costs.cel * sum(Rel) + costs.cbat * sum(Rbat); %Cost taking into account ageing of resources
    obj = obj + costs.cfc * (sum(c(fcindices).^2 + 100* c(fcindicesc).^2)) + costs.cel * sum(c(elindices).^2) + ...
        costs.cbat * (sum(c(batteryindices).^2 + 100*c(batteryindicesc).^2));
    obj = obj + sum((s(qsinds)).^2)*1e4;  
    %% Solve 
    
    %disp('Solving Optimization Problem (Subproblem)')
    options = sdpsettings('solver','gurobi','gurobi.Method',2,'verbose',1,'debug',1, 'gurobi.QCPDual', 1);%,'gurobi.BarQCPConvTol',1e-5);
    %options = sdpsettings('solver','mosek','verbose',0);
    sol = optimize(cons, obj,options);
    prob = sol.problem;
    
    if prob ~= 0
        disp('!!!!!!!!!!! Error Code !!!!!!!!!!!!')
        disp(prob)
        disp('Breakpoint')
        disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    end
    
    % Check Loadflow
    Errors = zeros(T,1); ErrorsS = zeros(T,1); ErrorsI = zeros(T,1);
    E_0 = ones(N,1);
    sol.ps = value(s(psinds)); sol.qs = value(s(qsinds)); sol.Pb = value(c(batteryindices));
    sol.Pfc = value(c(fcindices)); sol.Pel =  - value(c(elindices)); 
    sol.V = Vbar + value(Avsparse * (c - cbar)); sol.I = Ibar + value(Aisparse * (c - cbar)); sol.qc = value(c(qindices));
    sol.mh = value(mhcons - mhprod); sol.ph = value(ph);
    sol.Eb = value(Eb);

    for t=1:T
        Sinj = zeros(N,1);
        startposqc = (t-1) * (NC-3) + 1;
        startposdc = (t-1) * NC + 1;
        Sinj(1) = sol.ps(t) + 1i * sol.qs(t);
        Sinj(loads.indices(1)) = -loads.data{1}(t,scenario)*1/Sb;
        Sinj(loads.indices(2)) = -loads.data{2}(t,scenario)*1/Sb;
        Sinj(uncontrollable.indices(1)) = uncontrollable.data{1}(t,scenario)*1/Sb;
        Sinj(uncontrollable.indices(2)) = uncontrollable.data{2}(t,scenario)*1/Sb;
        Sinj(idxbat) = sol.Pb(t) + 1i * sol.qc(startposqc);
        Sinj(idxh) = sol.Pfc(t)-sol.Pel(t) + 1i * sol.qc(startposqc+1);
        startnode = (t-1) * N + 1;
        Vt = sol.V(startnode:startnode+N-1);
        [J,E,Snew,n_iter] = NR_polar_sol(Sinj,Vt,YY,E_0,idx,paramsNR);
        Sbars((t-1)*N+1:t*N) = Snew;
        Lbar((t-1)*2+1:2*t) = [sum(real(Snew));sum(imag(Snew))];
        Errors(t) = norm(abs(E)-Vt);
        psbar(t) = real(Snew(1)); qsbar(t) = imag(Snew(1));
        sbar(2*(t-1) + 1: 2*(t-1) +2) = [real(Snew(1)); imag(Snew(1))];
        ErrorsS(t) = abs(psbar(t)-sol.ps(t));
        
        Vbar(startnode:startnode+N-1) = E;
        cbar(startposdc:startposdc+NC-1) = [sol.Pb(t);sol.Pfc(t);-sol.Pel(t);sol.qc(startposqc:startposqc+1)];
        
        Itemp = zeros(NL,1);
        for i=1:NL
            idx1 = LP(i,1); idx2 = LP(i,2);
            Itemp(i) = YYL(idx1,idx2)*(E(idx1)-E(idx2)) + YYT(idx1,idx2)*E(idx1);
        end
        Ibar((t-1)*NL+1:t*NL) = abs(Itemp);
        ErrorsI(t) = norm(abs(Itemp) - abs(sol.I((t-1)*NL+1:t*NL)));
    end
   disp('LF:  '); disp(max(Errors)); disp(max(ErrorsS)); disp(max(ErrorsI))
    if max(Errors) <= tol 
        converged = 1;
    end
end

disp('LF:  '); disp(max(Errors));
disp('Nb of iterations:'); disp(niter);
sol.ps = psbar; sol.qs = qsbar; 
sol.V = abs(Vbar); sol.I = abs(Ibar); 
sol.pb = sol.Pb;
sol.dphslack = value(dphslack);

iterationcost = value(obj);



%% Plotting
plotting = 1;
if plotting ==1
    time = [linspace(0.25,24,24*4)];
    figure(1)
    for i =1:14
        subplot(4,4,i)
        Vs = reshape(sol.V,N,[]);
        plot(time,Vs(i,:))
        title('Nodal Voltages')
    end

    figure(2)
    ActiveInj = cell(3,1); ActiveInj{1} = 'Battery'; ActiveInj{2} = 'Fuel Cell'; ActiveInj{3} =  'Electrolyzer';
    subplot(3,3,1)
    pbs = sol.Pb;
    psum = pbs*Sb/1000;
    plot(time,pbs*Sb/1000)
    title(ActiveInj(1))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    subplot(3,3,2)
    hold on
    pfcs = sol.Pfc;
    psum=psum + pfcs*Sb/1000;
    plot(time,pfcs*Sb/1000)
    title(ActiveInj(2))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    subplot(3,3,3)
    hold on
    pels = sol.Pel;
    psum=psum - pels*Sb/1000;
    plot(time,pels*Sb/1000)
    title(ActiveInj(3))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    
    for i=1:length(loads.indices)
        subplot(3,3,3+i)
        L = loads.data{i};
        Lds = L(1:T,scenario);
        psum=psum - Lds/1000;
        plot(time,-Lds/1000)
        xlabel('Time (hours)')
        title('Load')
        ylabel('Power [kW]')
    end
    for i=1:length(uncontrollable.indices)
        subplot(3,3,5+i)
        U = uncontrollable.data{i};
        Us = U(1:T,scenario);
        psum=psum + Us/1000;
        plot(time,Us/1000)
        xlabel('Time (hours)')
        title('PV Generation ')
        ylabel('Power [kW]')
    end
    subplot(3,3,8)
    hold on
    pss = sol.ps;
    psum=psum + pss*Sb/1000;
    plot(time,pss*Sb/1000)
    hold off
    title('Slack Power')
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    
    subplot(3,3,9)
    Ebs = [sol.Eb(1:end-1)];
    plot(time,Ebs*Sb/1000)
    title('Battery Energy')
    xlabel('Time (hours)')
    ylabel('Energy [kWh]')
    
    figure(3)
    hold on
    phs = [sol.ph(1:end-1)];
    plot(time,phs*pb);
    title('Hydrogen Storage')
    xlabel('Time (hours)')
    ylabel('Pressure [bar]')
    
end

end
