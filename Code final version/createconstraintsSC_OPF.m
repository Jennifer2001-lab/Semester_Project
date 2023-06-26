function constraints = createconstraintsSC_OPF(topology, T, xbar, vbar, Ibar, Sbar, limits)
% CREATECONSTRAINTS SC OPF
% Return the constraints representing the Sensitivity Coefficient grid model for all
% timesteps 1..T 

%% Get the topology information
YY = topology.YY; YYL = topology.YYL; YYT = topology.YYL; indices = topology.indices; LP = topology.LP;
Ninj = size(indices.Active,2) + size(indices.Reactive,2); Nlines = size(YY,1)-1; Nnodes = size(YY,1); 
alphap = zeros(1,Nnodes); alphaq = zeros(1,Nnodes); % Voltage dependency of injections (not considered here)
Res_nodes_no_slack = [2:Nnodes]; slack=1; nph=1; % Parameters for SC computation

nineq = T*(2*(Nnodes-1) + 2*Nlines); neq = 2*T; % Don't check slack node, fixed to 1
nvar = (Ninj + 4)*T; %Controllable injections + slack (+/- for active power)

Aeq = zeros(neq, nvar); beq = zeros(neq,1); 
Aineq = zeros(nineq, nvar); bineq = zeros(nineq,1);

%% Compute the sensitivity coefficients
Av = zeros(Nlines*T,Ninj*T); Ai = Av; Ala = zeros(T,Ninj*T); Alr = Ala;
for t=1:T
    E = [1;vbar(t:T:(Nlines-1)*T+t)]; S = Sbar(t:T:(Nnodes-1)*T+t);
    % Voltage Sensitivity
    [K_p,K_com_p,K_q,K_comp_q]=Coeffs_Voltage_Alpha(YY,S,E,Res_nodes_no_slack,slack,nph,alphap,alphaq,ones(Nnodes,1));
    Kv = [K_p(2:end,indices.Active), K_q(2:end,indices.Reactive)];
    % Current Sensitivity
    [Icoeff, Icoeff_complex, Icoeffq,Icoeffq_complex, Currents]=Coeffs_Currents(YYL,YYT,E,K_com_p,K_comp_q,nph);
    % Write in easier format for constraints
    Aim = zeros(Nlines,(Nnodes-1)); % Sensitivities for one! sides of branches + wrt p and q
    for l=1:Nnodes-1
        Ic = cell2mat(Icoeff{l});
        Icq = cell2mat(Icoeffq{l});
        for i=1:Nlines
            idx1 = LP(i,1); idx2 = LP(i,2);
            Aim(i,l) = Ic(idx1,idx2);
           % Aim(i+NL,l) = Ic(idx2,idx1);
            Aim(i,l+Nnodes-1) = Icq(idx1,idx2);
           % Aim(i+NL,l+N-1) = Icq(idx2,idx1);
        end
    end
    Ki = [Aim(:,indices.Active), Aim(:,Nlines+indices.Reactive)];
    % Sensitivity of the losses
    [C_rp , C_rq, C_xp, C_xq] =Coeffs_Losses(YY, E, K_com_p, K_comp_q, Res_nodes_no_slack);
    Kla = [C_rp(:,indices.Active), C_rq(:,indices.Reactive)]; Klr =  [C_xp(:,indices.Active), C_xq(:,indices.Reactive)];
    Av(t:T:(Nlines-1)*T+t,t:T:(Ninj-1)*T+t) = Kv;
    Ai(t:T:(Nlines-1)*T+t,t:T:(Ninj-1)*T+t) = Ki;
    Ala(t,t:T:(Ninj-1)*T+t) = Kla; Alr(t,t:T:(Ninj-1)*T+t) = Klr;
end


%% Write Equations
startineq = 1; starteq = 1;
% Voltage limits
Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) = Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) + Av;
bineq(startineq:startineq + Nlines*T-1) = limits.vmax - abs(vbar) + Av * xbar(1:Ninj*T);
startineq = startineq + Nlines*T;
Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) = Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) - Av; bineq(startineq:startineq + Nlines*T-1) = -limits.vmin + abs(vbar) - Av * xbar(1:Ninj*T);
startineq = startineq + Nlines*T;

% Current limits
Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) = Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) + Ai; bineq(startineq:startineq + Nlines*T-1) = limits.Imax - Ibar + Ai * xbar(1:Ninj*T);
startineq = startineq + Nlines*T;
Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) = Aineq(startineq:startineq + Nlines*T-1, 1:Ninj*T) - Ai; bineq(startineq:startineq + Nlines*T-1) = -limits.Imin + Ibar - Ai * xbar(1:Ninj*T);
startineq = startineq + Nlines*T;

% Update Slack ACTIVE
Aeq(starteq:starteq+T-1, Ninj*T+1:Ninj*T+T) = Aeq(starteq:starteq+T-1, Ninj*T+1:Ninj*T+T) + eye(T); % Positive Active slack power
Aeq(starteq:starteq+T-1, Ninj*T+T+1:Ninj*T+2*T) = Aeq(starteq:starteq+T-1, Ninj*T+T+1:Ninj*T+2*T) - eye(T); % Negative Active slack power
Aeq(starteq:starteq+T-1, 1:Ninj*T) = Aeq(starteq:starteq+T-1, 1:Ninj*T) - Ala; % Losses sensitivity of injections
Aeq(starteq:starteq+T-1, 1:T) = Aeq(starteq:starteq+T-1, 1:T) + eye(T); % Compensation from change of injection of battery
Aeq(starteq:starteq+T-1, T+1:2*T) = Aeq(starteq:starteq+T-1, T+1:2*T) + eye(T); % Compensation from change of injection of fuel cell
beq(starteq:starteq+T-1) = xbar(Ninj*T+1:Ninj*T+T) - xbar(Ninj*T+T+1:Ninj*T+2*T); % Previous active slack power
beq(starteq:starteq+T-1) = beq(starteq:starteq+T-1) - Ala*xbar(1:Ninj*T);% Losses sensitivity of injections
beq(starteq:starteq+T-1) = beq(starteq:starteq+T-1) + xbar(1:T) + xbar(T+1:2*T); % Previous injections of active power
starteq = starteq + T;
% Update Slack REACTIVE
Aeq(starteq:starteq+T-1, Ninj*T+2*T+1:Ninj*T+3*T) = Aeq(starteq:starteq+T-1, Ninj*T+2*T+1:Ninj*T+3*T) + eye(T); % Positive Reactive slack power
Aeq(starteq:starteq+T-1, Ninj*T+3*T+1:Ninj*T+4*T) = Aeq(starteq:starteq+T-1, Ninj*T+3*T+1:Ninj*T+4*T) - eye(T); % Negative Reactive Slack Power
Aeq(starteq:starteq+T-1, 1:Ninj*T) = Aeq(starteq:starteq+T-1, 1:Ninj*T) - Alr; % Losses sensitivity of injections
Aeq(starteq:starteq+T-1, 2*T+1:3*T) = Aeq(starteq:starteq+T-1, 2*T+1:3*T) + eye(T); % Compensation from change of injection of battery
Aeq(starteq:starteq+T-1, 3*T+1:4*T) = Aeq(starteq:starteq+T-1, 3*T+1:4*T) + eye(T); % Compensation from change of injection of fuel cell
beq(starteq:starteq+T-1) = xbar(Ninj*T+2*T+1:Ninj*T+3*T) - xbar(Ninj*T+3*T+1:Ninj*T+4*T); % Previous reactive slack power
beq(starteq:starteq+T-1) = beq(starteq:starteq+T-1) - Alr*xbar(1:Ninj*T);
beq(starteq:starteq+T-1) = beq(starteq:starteq+T-1) + xbar(2*T+1:3*T) + xbar(3*T+1:4*T); % Previous injections of active power

constraints.Aeq = Aeq; constraints.beq = beq;
constraints.Aineq = Aineq; constraints.bineq = bineq;
constraints.Av = Av; constraints.Ai = Ai; constraints.Ala = Ala; constraints.Alr = Alr;

end



