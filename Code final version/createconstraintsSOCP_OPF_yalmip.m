function cons = createconstraintsSOCP_OPF_yalmip(topology, T, Vnodes, Plines, Qlines, flines, Pinj, Qinj)
% CREATECONSTRAINTSRELAXEDOPF
% Return the constraints representing the Relaxed OPF grid model for all
% timesteps 1..T 
% Return constraints in yalmip form, for one scenario and all considered
% timesteps
%% Get the topology information
R = topology.R; X = topology.X; B = topology.B;
G = topology.downstream_lines; Up = topology.upstream_nodes; Down = topology.downstream_nodes;
injections = topology.Inj; 
Nlines = size(R,1); Nnodes = Nlines + 1;

cons = [];
%% Write Equations
% Loop over lines to write constraints
cons = [cons, Vnodes(1:T) == ones(T,1)]; % Slack node voltage

for l=1:Nlines

    %%%% Active Power Flow %%%%
    % Active power at top of line l
    LHS = Plines((l-1)*T+1:l*T); RHS = zeros(size(LHS));
    % Active Power Injection at the end of line l
    LN = Down(l);
    if injections(LN) > 0
        c = injections(LN);
        RHS = RHS + Pinj((c-1)*T+1:c*T);
    end
    % Lines downstream from line l
    for m=1:Nlines
        if G(l,m) == 1
            RHS = RHS + Plines((m-1)*T+1:m*T);
        end
    end
    % Losses
    RHS = RHS + flines((l-1)*T+1:l*T) * R(l);
    cons = [cons, LHS == RHS];
    
    %%%% Reactive Power Flow %%%%
    % Reactive power at top of line l
    LHS = Qlines((l-1)*T+1:T*l); RHS = zeros(size(LHS));
    % Reactive Power Injection at the end of line l
    if injections(l) > 0
        c = injections(l);
        RHS  = RHS + Qinj((c-1)*T+1:c*T);
    end
    % Lines downstream from line l
    for m=1:Nlines
        if G(l,m) == 1
            RHS = RHS + Qlines((m-1)*T+1:m*T);
        end
    end
    % Shunt Injections
    n_up = Up(l);
    n_down = Down(l);
    RHS = RHS - Vnodes((n_down-1)*T+1:n_down*T) * B(l)/2;
    RHS = RHS - Vnodes((n_up-1)*T+1:n_up*T) * B(l)/2;
    % Losses
    RHS = RHS + flines((l-1)*T+1:l*T) * X(l);
    cons = [cons, LHS == RHS];

    %%%% Voltage %%%%
    % Nodal Voltage
    LHS = Vnodes((n_down-1) * T + 1:T * n_down); RHS = zeros(size(LHS));
    % Upstream Nodal Voltage
    RHS = RHS + Vnodes((n_up-1) * T + 1:T * n_up);
    % Active Drop
    RHS = RHS - 2*R(l)* Plines((l-1) * T+1:T * l);
    % Reactive Drop
    RHS = RHS - 2*X(l)* Qlines((l-1) * T+1:T * l);
    % Shunt Drop
    RHS = RHS - X(l) * B(l) * Vnodes((n_up-1) * T + 1:T * n_up);
    % Resistive Drop
    RHS = RHS + (X(l)^2 + R(l)^2) * flines((l-1) * T+1:T * l);
    
    cons = [cons, LHS == RHS];

    % Relaxed Ampacity Limit
    for t=1:T
        cons = [cons, cone([Vnodes((n_down-1)*T+t)+flines((l-1)*T+t);[2*Plines((l-1)*T+t);2*Qlines((l-1)*T+t);Vnodes((n_down-1)*T+t)-flines((l-1)*T+t)]])];
    end

end
end



