close all;
load MicroGrid.mat
load Resourceparameters.mat
load microloadsreal.mat

loads.data = Loads_Scenarios; loads.indices = load_indices;
uncontrollable.data = UncGen_Scenarios; uncontrollable.indices = uncon_indices;

state.SOCH = 0.5; state.SOCB = 0.5; % Initial battery level, initial hydrogen pressure
T = 24*4; S = 10; Sb = 600000; Vb = 400;

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
cbar = zeros(T*NC,S); Lbar = zeros(2*T,S); sbar = Lbar;
qsbar = zeros(T,S); Vbar = zeros(N*T,S); Ibar = zeros(NL*T,S); 
Sbars = zeros(T*N,S); psbar = zeros(T,S);

for t=1:T  
    Sinj = zeros(N,S);
    for scenario=1:S
        Sinj(loads.indices(1),scenario) = -loads.data{1}(t,scenario)*1/Sb;
        Sinj(loads.indices(2),scenario) = -loads.data{2}(t,scenario)*1/Sb;
        Sinj(uncontrollable.indices(1),scenario) = uncontrollable.data{1}(t,scenario)*1/Sb;
        Sinj(uncontrollable.indices(2),scenario) = uncontrollable.data{2}(t,scenario)*1/Sb;
    

        % Compute LoadFlow solution of initial state
        [J,E,Slf,n_iter] = NR_polar_sol(Sinj(:,scenario),E_star,YY,E_0,idx,paramsNR);
        Slf_(:,scenario)=Slf;
        E_(:,scenario)=E;
        Sbars((t-1)*N+1:t*N, scenario) = Slf_(:,scenario);
        Lbar((t-1)*2+1:2*t, scenario) = [sum(real(Slf_(:,scenario)));sum(imag(Slf_(:,scenario)))];
    
        pt(:,scenario) = real(Slf_(:,scenario)); 
        qt(:,scenario) = imag(Slf_(:,scenario)); 
        sbar(2*(t-1)+1:2*(t-1) + 2,scenario) = [pt(1,scenario);qt(1,scenario)];
        psbar(t,scenario) = pt(1,scenario); qsbar(t,scenario) = qt(1,scenario); 
    
        Vbar((t-1)*N+1:t*N,scenario) = E_(:,scenario); 
        Itemp = zeros(NL,S);
    
        for i=1:NL
            idx1 = LP(i,1); idx2 = LP(i,2);
            Itemp(i,scenario) = YYL(idx1,idx2)*(E(idx1)-E(idx2)) + YYT(idx1,idx2)*E(idx1);
        end
    end
    Ibar((t-1)*NL+1:t*NL,scenario) = abs(Itemp(:,scenario));
end

nmax = 10; % Limit the number of iterations
alphap = zeros(1,N); 
alphaq = zeros(1,N); % Voltage dependency of injections (not considered here)
Res_nodes_no_slack = [2:N]; slack=1; nph=1; % Parameters for SC computation, grid connection is considered as slack node



while converged == 0 && niter <= nmax
    niter = niter+1; 
    %% Compute Sensitivity Coefficients
    for t=1:T
        for scenario=1:S
            E(:,scenario) = Vbar((t-1)*N+1:t*N,scenario); 
            St(:,scenario)= Sbars((t-1)*N+1:t*N,scenario); 
        
              % Voltage Sensitivity
             [K_p,K_com_p,K_q,K_comp_q]=Coeffs_Voltage_Alpha(YY,St(:,scenario),E(:,scenario),Res_nodes_no_slack,slack,nph,alphap,alphaq,E_0);
             % Current Sensitivity
             [Icoeff, Icoeff_complex, Icoeffq,Icoeffq_complex, Currents]=Coeffs_Currents(YYL,YYT,E(:,scenario),K_com_p,K_comp_q,nph);
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
            [C_rp , C_rq, C_xp, C_xq] =Coeffs_Losses(YY, E(:,scenario), K_com_p, K_comp_q, Res_nodes_no_slack);
        
         % Save SC for every timestep and every scenario
             SC.Al{t,scenario} = [C_rp(idxbat-1),C_rp(idxh-1), C_rp(idxh-1), C_rq(idxbat-1),C_rq(idxh-1); ...
                 C_xp(idxbat-1), C_xp(idxh-1), C_xp(idxh-1), C_xq(idxbat-1), C_xq(idxh-1)];
             SC.Avm{t,scenario} = [K_p(:,idxbat-1),  K_p(:,idxh-1), K_p(:,idxh-1), K_q(:,idxbat-1), K_q(:,idxh-1)];
             SC.Aim{t,scenario} = [Aim(:,idxbat-1),  Aim(:,idxh-1), Aim(:,idxh-1), Aim(:,idxbat-1 + NL),  Aim(:,idxh-1 + NL)];
        end
    end 
    
     Vbar = abs(Vbar);
    
    %% Define subproblem variables
    Dslackpos = sdpvar(T,S); 
    Dslackneg = sdpvar(T,S);
    dispatch=sdpvar(T,1);
    mhprod = sdpvar(T,S); 
    mhcons = sdpvar(T,S);
    ph = sdpvar(T+1,S); 
    Eb = sdpvar(T+1,S); 
    c = sdpvar(NC*T,S); % Active and reactive together
    s = sdpvar(2*T,S); 
    dL = sdpvar(2*T,S);
    Rbat = sdpvar(T-1,S); Rfc = sdpvar(T-1,S); Rel = sdpvar(T-1,S); % Take into account ramping (absolute values)
    
    
     %% Define constraints
     cons = [];
      cons = [cons, Dslackpos >= 0, Dslackneg >= 0, Dslackneg <= 1, Dslackpos <= 1];
     % Impose grid constraints using SC
     
     for scenario=1:S
         % Voltage limits
          Avsparse = sparse(blkdiag(SC.Avm{:,scenario}));
          cons = [cons, Avsparse * (c(:,scenario) - cbar(:,scenario)) >= vmin - Vbar(:,scenario)]; % c-cbar represents [delta_p, delta_q]
          cons = [cons, Avsparse * (c(:,scenario) - cbar(:,scenario)) <= vmax - Vbar(:,scenario)];
     
         % Current Limits
         Aisparse = sparse(blkdiag(SC.Aim{:,scenario}));
         cons = [cons, Aisparse * (c(:,scenario) - cbar(:,scenario)) >= -Imax - Ibar(:,scenario)];
         cons = [cons, Aisparse * (c(:,scenario) - cbar(:,scenario)) <= Imax - Ibar(:,scenario)];
          % Compute losses based on injections
         Alsparse = sparse(blkdiag(SC.Al{:,scenario}));
          cons = [cons, dL(:,scenario) == Alsparse * (c(:,scenario) - cbar(:,scenario))];
     
          % Compute power at the slack node
          Acol = [1,1,1,0,0;0,0,0,1,1]; Acolcell = repmat({Acol},1,T);
          %electrolyser has no reactove power, p and q 
          Acolsparse = sparse(blkdiag(Acolcell{:})); %s is slack, it's has two columns of p and q, dL is the losses 

           cons = [cons, (s(:,scenario)-sbar(:,scenario)) == dL(:,scenario) - Acolsparse * (c(:,scenario)-cbar(:,scenario)), dL<=1, dL>= -1, s >= -1, s <= 1];
     
          % Slack node
          cons = [cons, s(psinds,scenario) <= resources.Pgmax, s(psinds,scenario) >= resources.Pgmin];
     
          % Storage constraints
           cons = [cons, ph(1,scenario)== state.SOCH * resources.phmax]; % Storage constraints 
           cons = [cons, Eb(1,scenario) == state.SOCB * resources.Ebmax];
          % Resource constraints
          cons = [cons, c(fcindices,scenario) <= resources.PFCmax, c(fcindices,scenario) >= resources.PFCmin];
          cons = [cons, c(elindices,scenario) >= -resources.PELmax, c(elindices,scenario) <= -resources.PELmin];
          cons = [cons, c(fcindicesc,scenario) <= resources.QFCmax, c(fcindicesc,scenario) >= resources.QFCmin]; 
          cons = [cons, c(batteryindices,scenario) <= resources.Pbmax, c(batteryindicesc,scenario) <= resources.Qbmax  ,...
                resources.Pbmin <= c(batteryindices,scenario), resources.Qbmin <= c(batteryindicesc,scenario)]; 
           cons = [cons, resources.kh_fc * c(fcindices,scenario) == mhcons(:,scenario), -resources.kh_el * c(elindices,scenario) == mhprod(:,scenario)];
     
          % Battery State of Charge
          cons = [cons, Eb(2:end,scenario) == Eb(1:end-1,scenario) - c(batteryindices(1:T),scenario) *1/4];% Battery energy level (assuming timesteps of 15 minutes)
           cons = [cons, Eb(:,scenario) >= resources.Ebmin, Eb(:,scenario) <= resources.Ebmax]; %Storage must remain above minimum level and below max
           % Hydrogen Storage
           cons = [cons, ph(2:end,scenario) == ph(1:end-1,scenario) - resources.Ktank * (mhcons(:,scenario) - mhprod(:,scenario))*1/4]; % Hydrogen tank pressure
           cons = [cons, ph(:,scenario) >= resources.phmin , ph(:,scenario) <= resources.phmax ]; % Hydrogen tank pressure limits
     
          cons = [cons, Rbat(:,scenario) >= 0, Rbat(:,scenario) >= c(batnext,scenario) - c(batprev,scenario),...
               Rbat(:,scenario) >= c(batprev,scenario) - c(batnext,scenario)];

          cons = [cons, Rfc(:,scenario) >= 0, Rfc(:,scenario) >= c(fcnext,scenario) - c(fcprev,scenario),...
               Rfc(:,scenario) >= c(fcprev,scenario) - c(fcnext,scenario)];
          cons = [cons, Rel(:,scenario)>=0, Rel(:,scenario) >= c(elnext,scenario) - c(elprev,scenario),...
               Rel(:,scenario) >= c(elprev,scenario) - c(elnext,scenario)];

           % Slack Power and Limit on Power Factor
            cons = [cons, s(psinds,scenario) == dispatch + Dslackpos(:,scenario) - Dslackneg(:,scenario)];
            %cons = [cons, abs(dispatch) >= s(qsinds,scenario) * tan(pi/2-resources.theta_max), abs(dispatch) >= - s(qsinds,scenario) * tan(pi/2-resources.theta_max)]; % Power Factor Constraint
     end
       
     
      %% Write Objective
    obj=0;
    %obj = obj + dphslack *  costs.cslack;
    %obj = obj +  1e5 * ((sum(Dslackpos.^2) + sum(Dslackneg.^2)));  % cost for deviating from dispatch %%  
%     obj = 0;
%     for t=1:T
%         obj = obj + 1e5 * sum(abs(Dslackpos(t,:)) + abs(Dslackneg(t,:)));
%     end
%     %obj = obj + costs.cfc * sum(Rfc) + costs.cel * sum(Rel) + costs.cbat * sum(Rbat); %Cost taking into account ageing of resources
%     obj = obj + costs.cfc * (sum(c(fcindices).^2 + 100* c(fcindicesc).^2)) + costs.cel * sum(c(elindices).^2) + ...
%         costs.cbat * (sum(c(batteryindices).^2 + 100*c(batteryindicesc).^2));
    %obj = obj + sum((s(qsinds)).^2)*1e4;  
    
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
    
    
    %% Check Loadflow
    Errors = zeros(T,S); ErrorsS = zeros(T,S); ErrorsI = zeros(T,S);
    E_0 = ones(N,1);
    sol.ps(:,scenario) = value(s(psinds,scenario)); sol.qs(:,scenario) = value(s(qsinds,scenario)); sol.Pb(:,scenario) = value(c(batteryindices,scenario));
    sol.Pfc(:,scenario) = value(c(fcindices,scenario)); sol.Pel(:,scenario) =  - value(c(elindices,scenario)); 
    sol.V = Vbar + value(Avsparse * (c - cbar)); sol.I = Ibar + value(Aisparse * (c - cbar)); sol.qc(:,scenario) = value(c(qindices,scenario));
    sol.mh = value(mhcons - mhprod); sol.ph = value(ph);
    sol.Eb = value(Eb);

    for t=1:T
        for scenario=1:S
            Sinj = zeros(N,S);
            startposqc = (t-1) * (NC-3) + 1;
            startposdc = (t-1) * NC + 1;
            Sinj(1) = sol.ps(t) + 1i * sol.qs(t);
            Sinj(loads.indices(1),scenario) = -loads.data{1}(t,scenario)*1/Sb;
            Sinj(loads.indices(2),scenario) = -loads.data{2}(t,scenario)*1/Sb;
            Sinj(uncontrollable.indices(1),scenario) = uncontrollable.data{1}(t,scenario)*1/Sb;
            Sinj(uncontrollable.indices(2),scenario) = uncontrollable.data{2}(t,scenario)*1/Sb;
            Sinj(idxbat,scenario) = sol.Pb(t) + 1i * sol.qc(startposqc);
            Sinj(idxh,scenario) = sol.Pfc(t)-sol.Pel(t) + 1i * sol.qc(startposqc+1);
            startnode = (t-1) * N + 1;
            Vt = sol.V(startnode:startnode+N-1,scenario);
            [J,E,Snew,n_iter] = NR_polar_sol(Sinj(:,scenario),Vt,YY,E_0,idx,paramsNR);
            Sbars((t-1)*N+1:t*N,scenario) = Snew;
            Lbar((t-1)*2+1:2*t,scenario) = [sum(real(Snew));sum(imag(Snew))];
            Errors(t,scenario) = norm(abs(E)-Vt);
            psbar(t,scenario) = real(Snew(1)); qsbar(t,scenario) = imag(Snew(1));
            sbar(2*(t-1) + 1: 2*(t-1) +2,scenario) = [real(Snew(1)); imag(Snew(1))];
            ErrorsS(t,scenario) = abs(psbar(t,scenario)-sol.ps(t));
        
            Vbar(startnode:startnode+N-1,scenario) = E;
            cbar(startposdc:startposdc+NC-1,scenario) = [sol.Pb(t,scenario);sol.Pfc(t,scenario);-sol.Pel(t,scenario);sol.qc(startposqc:startposqc+1,scenario)];
            Itemp = zeros(NL,S);
            for i=1:NL
                idx1 = LP(i,1); idx2 = LP(i,2);
                Itemp(i,scenario) = YYL(idx1,idx2)*(E(idx1)-E(idx2)) + YYT(idx1,idx2)*E(idx1);
            end
            Ibar((t-1)*NL+1:t*NL,scenario) = abs(Itemp(:,scenario));
            ErrorsI(t,scenario) = norm(abs(Itemp(:,scenario)) - abs(sol.I((t-1)*NL+1:t*NL,scenario)));
        end
    end
   disp('LF:  '); disp(max(Errors)); disp(max(ErrorsS)); disp(max(ErrorsI))
    if max(Errors) <= tol 
        converged = 1;
    end
end

%% Displaying errors
disp('LF:  '); disp(max(Errors));
disp('Nb of iterations:'); disp(niter);
sol.ps = psbar; sol.qs = qsbar; 
sol.V = abs(Vbar); sol.I = abs(Ibar); 
sol.pb = sol.Pb;

iterationcost = value(obj);

%% Plotting
scenario=1
plotting = 1;
if plotting ==1
    time = [linspace(0.25,24,24*4)];
    figure(1)
    for i =1:14
        subplot(4,4,i)
        Vs = reshape(sol.V(:,scenario),N,[]);
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
    pss = sol.ps(:,scenario);
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




