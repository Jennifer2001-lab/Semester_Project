function [objective, sol, dual_solution_sub] = solve_subproblem_sequence_matrix(loads, uncontrollable, grid, resources, costs, el_price, master_sol, T, S, pv_num, subprob_number, sequence_number, nSP)
    %% Unpack the Grid Data
    state.SOCH = 0.5; state.SOCB = 0.5; % Initial battery level, initial hydrogen pressure

    bv =grid.basevalues; % Base values
    Vb=bv(1); Sb=bv(2); pb=bv(5); Ib = Sb/(sqrt(3)*Vb);
    LP = grid.lineparameters; % Line parameters, indices, admittances, ampacities ...
    NL = grid.lines; N = grid.nodes; NC = 5; % Number of lines and nodes, controllable: battery (p and q), fuel cell (p and q), electrolyzer (q)
    idx = grid.node_indices; % nodal indices
    adm = grid.Admittance; % Admittance matrices
    YY = adm{1}; YYL=adm{2}; YYT=adm{3}; Ib=adm{4}; 
    idxbat = 5; idbat =1;
    idxh = 13; idfc = 2; idel = 3;
    Imax = repmat((LP(:,7)/Ib),T,1); vmin = 0.95; vmax = 1.05;
    limits.Imax = Imax; limits.Imin = -Imax; limits.vmax = vmax*ones(NL*T,1); limits.vmin = vmin*ones(NL*T,1);
    topology.Ninj = 4;
    topology.LP = LP;
    topology.indices.Active = [idxh-1,idxbat-1]; topology.indices.Reactive = topology.indices.Active;
    topology.YY = YY; topology.YYT = YYT; topology.YYL = YYL;
    paramsNR.tol=1e-8; paramsNR.n_max=1000;  % maximum number of iterations, convergence tolerance

    [slope_B, constant_B] = Linearize_quadratic(1, 6);
    [slope_FC, constant_FC] = Linearize_quadratic(1, 6);

    %% Solve initial LF
    Injbar = zeros(4*T,S); Ibar = zeros(NL*T,S); Sbar = zeros(N*T,S); Vbar = zeros((N-1)*T,S);
    pspbar = zeros(T,S); psnbar = zeros(T,S); qspbar = zeros(T,S); qsnbar = zeros(T,S);
    for s=1:S
        for t=1:T
            % Get the starting values for the controllable injections
            Sinj = zeros(N,1);
            Sinj(loads.indices(1)) = -10e3*loads.data(t,s)*1/Sb;
            Sinj(loads.indices(2)) = -5e3*loads.data(t,s)*1/Sb;
            Sinj(uncontrollable.indices(1)) = 10e3*uncontrollable.data{pv_num}(t,s)*1/Sb;
            Sinj(uncontrollable.indices(2)) = 10e3*uncontrollable.data{pv_num}(t,s)*1/Sb;
    
            % Compute LoadFlow solution of initial state
            [J,E,Slf,n_iter] = NR_polar_sol(Sinj,ones(N,1),YY,ones(N,1),idx,paramsNR);
            Sbar(t:T:(N-1)*T+t,s) = Slf;
            pspbar(t,s) = max(0, real(Slf(1))); psnbar(t,s) = -min(0, real(Slf(1)));
            qspbar(t,s) = max(0, imag(Slf(1))); qsnbar(t,s) = -min(0, imag(Slf(1)));
            Vbar(t:T:(N-2)*T+t,s) = E(2:end); 
            Itemp = zeros(NL,1);
            for i=1:NL
                idx1 = LP(i,1); idx2 = LP(i,2);
                Itemp(i) = YYL(idx1,idx2)*(E(idx1)-E(idx2)) + YYT(idx1,idx2)*E(idx1);
               % Itemp(i+NL) = YYL(idx2,idx1)*(E(idx2)-E(idx1)) +
               % YYT(idx2,idx1)*E(idx2); Ideally also check second side
            end
            Ibar(t:T:(NL-1)*T+t,s) = abs(Itemp);
        end
    end
    xbar = [Injbar;pspbar;psnbar;qspbar; qsnbar];
    converged = 0; niter = 0; nmax = 5; tol = 1e-6; 
    while converged == 0 && niter < nmax
        niter = niter+1; 
        %% Define the resource constraints
        ncomvar = 9 + 2*T + 2;
        nrecvar = 2 * T + 2 * (T+1) + 3 * (T-1) + 2*T + 2*T;
        ngridvar = topology.Ninj * T + 4*T;

        ncommoneq = 10;
        nreceq = (T+1) * 2 + T + 2*T;
        nrecineq = 2*(T+1) + 3 * 2 * (T-1) + 12*T + 6*T;
        % X is defined as:
        % [Ebmax, Sbmax, Ehmax, PFCmax, PELmax, dEh, dEhp, dEhm, eps, Dplus, Dminus, ufc, uel, Dsplus, Dsminus, Ebat, Ehydrogen, Ramping, Pbp, Pbn, Pfc, Pel, Injections, Ps+, Ps-, Qs+, Qs-]
        offsetEbmax = 1; offsetSbmax = 2; offsetEhmax = 3; offsetPFCmax = 4; offsetPELmax = 5; offsetdEh = 6; 
        offsetdEhp = offsetdEh + 1; offsetdEhm = offsetdEhp+1; offseteps = 7;
        offsetDplus = offseteps + 1; offsetDminus = offsetDplus + T; offsetfcbin = offsetDminus + T; offsetfelbin = offsetfcbin + 1;
        offsetCommon = offsetfelbin + 1;
        
        Acomeq = zeros(ncommoneq, ncomvar + S*(nrecvar + ngridvar)); bcomeq = zeros(ncommoneq,1);
        Aeq = []; Aineq = []; beq = []; bineq = []; ub = []; lb = []; f = [];
        %% First write the equations only concerning the common variables (for every scenario)
        neq=1;
        Acomeq(neq,offsetEbmax) = 1; bcomeq(neq) = master_sol.Ebmax;
        neq = neq + 1;
        Acomeq(neq,offsetSbmax) = 1; bcomeq(neq) = master_sol.Sbmax;
        neq = neq + 1;
        Acomeq(neq,offsetEhmax) = 1; bcomeq(neq) = master_sol.Ehmax;
        neq = neq + 1;
        Acomeq(neq,offsetPFCmax) = 1; bcomeq(neq) = master_sol.PFCmax;
        neq = neq + 1;
        Acomeq(neq,offsetPELmax) = 1; bcomeq(neq) = master_sol.PELmax;
        neq = neq + 1;
        Acomeq(neq,offsetdEh) = 1; bcomeq(neq) = master_sol.dEh(subprob_number);
        neq = neq + 1;
        Acomeq(neq,offsetdEhp) = 1; bcomeq(neq) = master_sol.dEhp(subprob_number);
        neq = neq + 1;
        Acomeq(neq,offsetdEhm) = 1; bcomeq(neq) = master_sol.dEhm(subprob_number);
        neq = neq + 1;
        Acomeq(neq,offsetfcbin) = 1; bcomeq(neq) = master_sol.fc_bin(:,subprob_number);
        neq = neq + 1;
        Acomeq(neq,offsetfelbin) = 1; bcomeq(neq) = master_sol.el_bin(:,subprob_number);
        neq = neq + 1;

        % Add the equations to the full system
        Aeq = [Aeq; Acomeq]; beq = [beq; bcomeq];
        
        ub = [ones(8,1); 10; resources.Pgmax*ones(2*T,1); ones(2,1)];
        lb = [zeros(9,1); zeros(2*T,1); zeros(2,1)];
        f = [zeros(8,1); 1e6; el_price(1:T)'; -1/3*el_price(1:T)'; zeros(2,1)];
        
        %% Loop over the scenarios and write the constraints specific to a scenario (eventually linked with the common variables)
        for s=1:S
            Areq = spalloc(nreceq, ncomvar + S*(nrecvar + ngridvar), nreceq*10); breq = zeros(nreceq,1);
            Arineq = spalloc(nrecineq, ncomvar + S*(nrecvar + ngridvar), nrecineq*10); brineq = zeros(nrecineq,1);
            % Offsets for variable and equation locations
            neq = 1;
            offsetDsplus = (s-1)*(ngridvar + nrecvar) + ncomvar + 1;
            offsetDsminus = offsetDsplus + T; offsetEbat = offsetDsminus + T; offsetEhydrogen = offsetEbat + T+1; 
            offsetRamping = offsetEhydrogen + T+1; offsetPbp = offsetRamping + 3*(T-1); offsetPbn = offsetPbp + T;
            offsetPfc = offsetPbn + T; offsetPel = offsetPfc + T;
            offsetInj = offsetPel + T; offsetPs_plus = offsetInj + 4*T; offsetPs_minus = offsetPs_plus + T; 
            offsetQsplus = offsetPs_minus + T; offsetQsminus = offsetQsplus + T;
        
            % Dispatch Constraint (Here Dsplus and Dsminus are the dispatch errors)
            Areq(neq:neq+ T-1, offsetDsplus:offsetDsminus-1) = Areq(neq:neq+ T-1, offsetDsplus:offsetDsminus-1) + eye(T);
            Areq(neq:neq+ T-1, offsetDsminus:offsetEbat-1) = Areq(neq:neq+ T-1, offsetDsminus:offsetEbat-1) - eye(T);
            Areq(neq:neq+ T-1, offsetPs_plus:offsetPs_plus+T-1) = Areq(neq:neq+ T-1, offsetPs_plus:offsetPs_plus+T-1) + eye(T);
            Areq(neq:neq+ T-1, offsetPs_minus:offsetPs_minus+T-1) = Areq(neq:neq+ T-1, offsetPs_minus:offsetPs_minus+T-1) - eye(T);
            % Set equal to common dispatch
            Areq(neq:neq+ T-1, offsetDplus:offsetDplus+T-1) = Areq(neq:neq+ T-1, offsetDplus:offsetDplus+T-1) - eye(T);
            Areq(neq:neq+ T-1, offsetDminus:offsetDminus+T-1) = Areq(neq:neq+ T-1, offsetDminus:offsetDminus +T-1) + eye(T);
            % Update the equation counter
            neq = neq + T;
            
            % Energy Storage Initializations
            % Initial states for energy storages (bat and hydr)
            Areq(neq,offsetEbat) = 1; breq(neq) = 0.5*resources.Ebmax;
            Areq(neq+1,offsetEhydrogen) = 1; breq(neq+1) = 0;
            neq = neq + 2;
            
            % Energy Storage Updates
            Areq(neq:neq + T-1, offsetEbat+1: offsetEbat+T) = Areq(neq:neq + T-1, offsetEbat+1: offsetEbat+T) + eye(T);
            Areq(neq:neq + T-1, offsetEbat: offsetEbat+T-1) = Areq(neq:neq + T-1, offsetEbat: offsetEbat+T-1) - eye(T);
            Areq(neq:neq + T-1, offsetPbp: offsetPbp+T-1) = Areq(neq:neq + T-1, offsetPbp: offsetPbp+T-1) - eye(T)*1/0.95;
            Areq(neq:neq + T-1, offsetPbn: offsetPbn+T-1) = Areq(neq:neq + T-1, offsetPbn: offsetPbn+T-1) + eye(T)*0.95;
            neq = neq + T;
            
            % Energy update for hydrogen (day, top)
            %Day
            Areq(neq:neq + T-1, offsetEhydrogen+1:offsetEhydrogen+T) = Areq(neq:neq + T-1, offsetEhydrogen+1:offsetEhydrogen+T) + eye(T);
            Areq(neq:neq + T-1, offsetEhydrogen:offsetEhydrogen+T-1) = Areq(neq:neq + T-1, offsetEhydrogen:offsetEhydrogen+T-1) - eye(T);
            Areq(neq:neq + T-1, offsetPfc:offsetPfc+T-1) = Areq(neq:neq + T-1, offsetPfc:offsetPfc+T-1) + eye(T) * 1/0.75;
            Areq(neq:neq + T-1, offsetPel:offsetPel+T-1) = Areq(neq:neq + T-1, offsetPel:offsetPel+T-1) - eye(T) * 0.75;
            neq = neq+T;
            
            % Power Balance for battery 
            offsetPbat = offsetInj;
            Areq(neq:neq + T-1, offsetPbp:offsetPbp+T-1) = Areq(neq:neq + T-1, offsetPbp:offsetPbp+T-1) + eye(T);
            Areq(neq:neq + T-1, offsetPbn:offsetPbn+T-1) = Areq(neq:neq + T-1, offsetPbn:offsetPbn+T-1) - eye(T);
            Areq(neq:neq + T-1, offsetPbat:offsetPbat+T-1) = Areq(neq:neq + T-1, offsetPbat:offsetPbat+T-1) - eye(T);
            neq = neq + T;
            % Power Balance at Hydrogen Node
            offsetPh2 = offsetPbat + T;
            Areq(neq:neq + T-1, offsetPfc:offsetPfc+T-1) = Areq(neq:neq + T-1, offsetPfc:offsetPfc+T-1) + eye(T);
            Areq(neq:neq + T-1, offsetPel:offsetPel+T-1) = Areq(neq:neq + T-1, offsetPel:offsetPel+T-1) - eye(T);
            Areq(neq:neq + T-1, offsetPh2:offsetPh2+T-1) = Areq(neq:neq + T-1, offsetPh2:offsetPh2+T-1) - eye(T);
            neq = neq + T;
            
            nineq=1;
            % Hydrogen Storage Limits
            
            Arineq(nineq:nineq+T-1,offsetEhydrogen:offsetEhydrogen+T-1) = Arineq(nineq:nineq+T-1,offsetEhydrogen:offsetEhydrogen+T-1) - eye(T);
            Arineq(nineq:nineq+T-1, offseteps) = -ones(T,1);
            Arineq(nineq:nineq+T-1, offsetdEhm) = ones(T,1);
            nineq=nineq+T;
            Arineq(nineq:nineq+T-1,offsetEhydrogen:offsetEhydrogen+T-1) = Arineq(nineq:nineq+T-1,offsetEhydrogen:offsetEhydrogen+T-1) + eye(T);
            Arineq(nineq:nineq+T-1, offseteps) = ones(T,1);
            Arineq(nineq:nineq+T-1, offsetdEhp) = ones(T,1);
            nineq = nineq+T;

            Arineq(nineq,offsetEhydrogen+T) = 1; Arineq(nineq,offsetdEh) = -1; Arineq(nineq,offseteps) = -1;
            nineq = nineq+1;
            Arineq(nineq,offsetEhydrogen+T) = -1; Arineq(nineq,offsetdEh) = 1; Arineq(nineq,offseteps) = -1;
            nineq = nineq+1;
            % Power Factor Constraint: skip for now
            % Resource Capability constraints
            % Unit commitment FC
            Arineq(nineq:nineq+T-1,offsetfcbin) = -ones(T,1) * 10;
            Arineq(nineq:nineq+T-1,offsetPfc:offsetPfc+T-1) = Arineq(nineq:nineq+T-1,offsetPfc:offsetPfc+T-1) + eye(T);
            nineq = nineq+T;
            Arineq(nineq:nineq+T-1,offsetPFCmax) = Arineq(nineq:nineq+T-1,offsetPFCmax) - ones(T,1);
            Arineq(nineq:nineq+T-1,offsetPfc:offsetPfc+T-1) = Arineq(nineq:nineq+T-1,offsetPfc:offsetPfc+T-1) + eye(T);
            nineq = nineq + T;
            Arineq(nineq:nineq+T-1,offsetPFCmax) = 0.1*ones(T,1);
            Arineq(nineq:nineq+T-1,offsetPfc:offsetPfc+T-1) = Arineq(nineq:nineq+T-1,offsetPfc:offsetPfc+T-1) - eye(T);
            Arineq(nineq:nineq+T-1,offsetfcbin) = ones(T,1) * 10;
            brineq(nineq:nineq+T-1) = 10*ones(T,1);
            nineq = nineq+T;
            % Unit commitment EL
            Arineq(nineq:nineq+T-1,offsetfelbin) = -ones(T,1) * 10;
            Arineq(nineq:nineq+T-1,offsetPel:offsetPel+T-1) = Arineq(nineq:nineq+T-1,offsetPel:offsetPel+T-1) + eye(T);
            nineq = nineq+T;
            Arineq(nineq:nineq+T-1,offsetPELmax) = Arineq(nineq:nineq+T-1,offsetPELmax) - ones(T,1);
            Arineq(nineq:nineq+T-1,offsetPel:offsetPel+T-1) = Arineq(nineq:nineq+T-1,offsetPel:offsetPel+T-1) + eye(T);
            nineq = nineq + T;
            Arineq(nineq:nineq+T-1,offsetPELmax) = 0.1*ones(T,1);
            Arineq(nineq:nineq+T-1,offsetPel:offsetPel+T-1) = Arineq(nineq:nineq+T-1,offsetPel:offsetPel+T-1) - eye(T);
            Arineq(nineq:nineq+T-1,offsetfelbin) = ones(T,1) * 10;
            brineq(nineq:nineq+T-1) = 10*ones(T,1);
            nineq = nineq+T;
        
            % PQ relationship
            offsetQbat = offsetPbat+2*T;
            for sl = 1:numel(slope_B)/2
                Arineq(nineq+(sl-1)*2*T:nineq+T-1+(sl-1)*2*T, offsetPbat:offsetPbat+T-1) = Arineq(nineq+(sl-1)*2*T:nineq+T-1+(sl-1)*2*T, offsetPbat:offsetPbat+T-1) - eye(T);
                Arineq(nineq+(sl-1)*2*T:nineq+T-1+(sl-1)*2*T, offsetQbat:offsetQbat+T-1) = Arineq(nineq+(sl-1)*2*T:nineq+T-1+(sl-1)*2*T, offsetQbat:offsetQbat+T-1) - slope_B(sl) * eye(T);
                Arineq(nineq+(sl-1)*2*T:nineq+T-1+(sl-1)*2*T, offsetSbmax) = -constant_B(sl);
                Arineq(nineq+(sl-1)*2*T+T:nineq+T-1+(sl-1)*2*T+T, offsetPbat:offsetPbat+T-1) = Arineq(nineq+(sl-1)*2*T+T:nineq+T-1+(sl-1)*2*T+T, offsetPbat:offsetPbat+T-1) + eye(T);
                Arineq(nineq+(sl-1)*2*T+T:nineq+T-1+(sl-1)*2*T+T, offsetQbat:offsetQbat+T-1) = Arineq(nineq+(sl-1)*2*T+T:nineq+T-1+(sl-1)*2*T+T, offsetQbat:offsetQbat+T-1) + slope_B(sl) * eye(T);
                Arineq(nineq+(sl-1)*2*T+T:nineq+T-1+(sl-1)*2*T+T, offsetSbmax) = -constant_B(sl);
            end
            nineq = nineq + sl*2*T;

            % Ramping
            % Battery
            Arineq(nineq:nineq+T-2,offsetRamping:offsetRamping+T-2) = Arineq(nineq:nineq+T-2,offsetRamping:offsetRamping+T-2) - eye(T-1);
            Arineq(nineq:nineq+T-2,offsetPbat+1:offsetPbat+T-1) = Arineq(nineq:nineq+T-2,offsetPbat+1:offsetPbat+T-1) - eye(T-1);
            Arineq(nineq:nineq+T-2,offsetPbat:offsetPbat+T-2) = Arineq(nineq:nineq+T-2,offsetPbat:offsetPbat+T-2) + eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetRamping:offsetRamping+T-2) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetRamping:offsetRamping+T-2) - eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetPbat+1:offsetPbat+T-1) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetPbat+1:offsetPbat+T-1) + eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetPbat:offsetPbat+T-2) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetPbat:offsetPbat+T-2) - eye(T-1);
            nineq = nineq + 2*(T-1);
            % Fuel Cell
            Arineq(nineq:nineq+T-2,offsetRamping+T-1:offsetRamping+T-1+T-2) = Arineq(nineq:nineq+T-2,offsetRamping+T-1:offsetRamping+T-1+T-2) - eye(T-1);
            Arineq(nineq:nineq+T-2,offsetPfc+1:offsetPfc+T-1) = Arineq(nineq:nineq+T-2,offsetPfc+1:offsetPfc+T-1) - eye(T-1);
            Arineq(nineq:nineq+T-2,offsetPfc:offsetPfc+T-2) = Arineq(nineq:nineq+T-2,offsetPfc:offsetPfc+T-2) + eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetRamping+T-1:offsetRamping+T-1+T-2) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetRamping+T-1:offsetRamping+T-1+T-2) - eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetPfc+1:offsetPfc+T-1) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetPfc+1:offsetPfc+T-1) + eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetPfc:offsetPfc+T-2) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetPfc:offsetPfc+T-2) - eye(T-1);
            nineq = nineq + 2*(T-1);
            % Electrolyzer
            Arineq(nineq:nineq+T-2,offsetRamping+2*(T-1):offsetRamping+2*(T-1)+T-2) = Arineq(nineq:nineq+T-2,offsetRamping+2*(T-1):offsetRamping+2*(T-1)+T-2) - eye(T-1);
            Arineq(nineq:nineq+T-2,offsetPel+1:offsetPel+T-1) = Arineq(nineq:nineq+T-2,offsetPel+1:offsetPel+T-1) - eye(T-1);
            Arineq(nineq:nineq+T-2,offsetPel:offsetPel+T-2) = Arineq(nineq:nineq+T-2,offsetPel:offsetPel+T-2) + eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetRamping+2*(T-1):offsetRamping+2*(T-1)+T-2) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetRamping+2*(T-1):offsetRamping+2*(T-1)+T-2) - eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetPel+1:offsetPel+T-1) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetPel+1:offsetPel+T-1) + eye(T-1);
            Arineq(nineq+T-1:nineq+T-1+T-2,offsetPel:offsetPel+T-2) = Arineq(nineq+T-1:nineq+T-1+T-2,offsetPel:offsetPel+T-2) - eye(T-1);
            nineq = nineq + 2*(T-1);
            
            %% Get the grid constraints (SC)
            grid_constraints = createconstraintsSC_OPF(topology, T, xbar(:,s), Vbar(:,s), Ibar(:,s), Sbar(:,s), limits);
            Agrideq = spalloc(size(grid_constraints.Aeq,1), ncomvar + S * (nrecvar + ngridvar), ngridvar*10);
            Agrideq(1:end,ncomvar + (s-1) * (nrecvar + ngridvar) + nrecvar+1:ncomvar+s * (nrecvar + ngridvar)) = grid_constraints.Aeq;
            Agridineq = spalloc(size(grid_constraints.Aineq,1), ncomvar + S * (nrecvar + ngridvar), ngridvar*10);
            Agridineq(1:end,ncomvar + (s-1) * (nrecvar + ngridvar) + nrecvar+1:ncomvar+s * (nrecvar + ngridvar)) = grid_constraints.Aineq;
        
            %% Also save the sensitivity coefficients
            AVS{s} = grid_constraints.Av; AIS{s} = grid_constraints.Ai; ALAS{s} = grid_constraints.Ala; ALRS{s} = grid_constraints.Alr;
        
            %% Add the constraints to the full system
            Aeq = sparse([Aeq;Areq;Agrideq]); beq = [beq;breq;grid_constraints.beq];
            Aineq = sparse([Aineq; Arineq; Agridineq]); bineq = [bineq; brineq;grid_constraints.bineq];
    
            %% Add bounds
            % [Ebmax, Sbmax, Ehmax, PFCmax, PELmax, dEh, eps, Dplus, Dminus, ufc, uel, Dsplus, Dsminus, Ebat, 
            % Ehydrogen, Ramping, Pbp, Pbn, Pfc, Pel, Injections, Ps+, Ps-, Qs+, Qs-]
            
            ub = [ub;resources.Pgmax*ones(2*T,1); resources.Ebmax*ones(T+1,1); ones(T+1,1); ones(3*(T-1),1); resources.Pbmax*ones(2*T,1);
                resources.PFCmax * ones(T,1); resources.PELmax*ones(T,1); resources.Pbmax*ones(T,1); resources.PFCmax*ones(T,1);zeros(2*T,1);
                resources.Pgmax*ones(4*T,1)];
            lb = [lb;zeros(2*T,1); resources.Ebmin*ones(T+1,1); -1*ones(T+1,1); zeros(3*(T-1),1); zeros(2*T,1);
                zeros(2*T,1); resources.Pbmin*ones(T,1); -resources.PELmax*ones(T,1);zeros(2*T,1);
                zeros(4*T,1)];
        
        
            %% Add Objective
            fscen = 1/S*[10 * max(el_price) * ones(2*T,1); zeros(2*T+2,1); costs.crbat*ones(T-1,1); costs.crfc*ones(T-1,1); costs.crel * ones(T-1,1);
               zeros(2*T,1); zeros(2*T,1); zeros(4*T,1); zeros(2*T,1);max(el_price)*ones(2*T,1)];
            f = [f;fscen];
        
        end
        f = 10*365/nSP *f;

        %% Create Structures to pass the model to gurobi
        vtype = repelem(char('C'),size(f,1))';
        sense = [repelem(char("="),size(Aeq,1))'; repelem(char("<"), size(Aineq,1))'];
        A = sparse([Aeq;Aineq]); b = [beq;bineq]; 
        %% Give the model to gurobi
        
        model.A = A; model.rhs = b; model.sense = sense; model.obj = f; model.lb = lb; model.ub = ub; 
        model.vtype = vtype;
        gurobi_write(model,'subproblem.lp')
        %% Solve the problem
        results = gurobi(model);
        x = results.x;
        duals = results.pi;


        dual_solution_sub.Ebmax = -duals(offsetEbmax);
        dual_solution_sub.Sbmax = -duals(offsetSbmax);
        dual_solution_sub.Ehmax = -duals(offsetEhmax);
        dual_solution_sub.PFCmax = -duals(offsetPFCmax);
        dual_solution_sub.PELmax = -duals(offsetPELmax);
        dual_solution_sub.dEh = -duals(offsetdEh);
        dual_solution_sub.dEhp = -duals(offsetdEhp);
        dual_solution_sub.dEhm = -duals(offsetdEhm);
        dual_solution_sub.fc_bin = -duals(offsetfcbin);
        dual_solution_sub.el_bin = -duals(offsetfelbin);


        status = results.status;

        % Check if optimization was successful
        if status =='OPTIMAL' 
            dual_solution_sub.warning = 0;
            objective = results.objval;
            sol.Ebmax=x(offsetEbmax);
            sol.Sbmax=x(offsetSbmax);
            sol.Ehmax=x(offsetEhmax);
            sol.PFCmax=x(offsetPFCmax);
            sol.PELmax=x(offsetPELmax);
            sol.dEh =x(offsetdEh);
            sol.dEhp =x(offsetdEhp);
            sol.dEhm =x(offsetdEhm);
            sol.el_bin=x(offsetfelbin);
            sol.fc_bin=x(offsetfcbin);

            sol.eps = x(offseteps);
            sol.dispatch = x(offsetDplus) - x(offsetDminus);

        else
            dual_solution_sub.warning = 1;
            % Infeasible or unbounded
            objective = 0;

            sol = master_sol;

            disp('Error in the optimization');

        end
        


        %% Check Loadflow
        Errors = zeros(T,S); ErrorsS = zeros(T,S); 
        E_0 = ones(N,1);
        sol.ps = zeros(T,S); sol.qs = sol.ps; sol.Pb = zeros(T,S); sol.Qb = zeros(T,S);
        sol.Pfc = zeros(T,S); sol.Pel = zeros(T,S); sol.Qfc = zeros(T,S);
        sol.Eb = zeros(T+1,S); sol.Eh= zeros(T+1,S);
        sol.V = zeros((N-1)*T,S); SLF = zeros(T,S); Losses = zeros(T,S);
        for s=1:S
            offsetDsplus = (s-1)*(ngridvar + nrecvar) + ncomvar + 1;
            offsetDsminus = offsetDsplus + T; offsetEbat = offsetDsminus + T; offsetEhydrogen = offsetEbat + T+1; 
            offsetRamping = offsetEhydrogen + T+1; offsetPbp = offsetRamping + 3*(T-1); offsetPbn = offsetPbp + T;
            offsetPfc = offsetPbn + T; offsetPel = offsetPfc + T;
            offsetInj = offsetPel + T; offsetPs_plus = offsetInj + 4*T; offsetPs_minus = offsetPs_plus + T; 
            offsetQsplus = offsetPs_minus + T; offsetQsminus = offsetQsplus + T;
            offsetPbat = offsetInj; offsetQbat = offsetPbat + 2*T; offsetQfc = offsetQbat+T;

            sol.Pb(:,s) = x(offsetPbat:offsetPbat+T-1); sol.Qb(:,s) = x(offsetQbat:offsetQbat+T-1); 
            sol.Pfc(:,s) = x(offsetPfc:offsetPfc+T-1); sol.Pel(:,s) = x(offsetPel:offsetPel+T-1); sol.Qfc(:,s) = x(offsetQfc:offsetQfc+T-1); 
            sol.Eb(:,s) = x(offsetEbat:offsetEbat+T); sol.Eh(:,s) = x(offsetEhydrogen:offsetEhydrogen+T); 
            sol.ps(:,s) = x(offsetPs_plus:offsetPs_plus+T-1) - x(offsetPs_minus:offsetPs_minus+T-1); sol.qs(:,s) = x(offsetQsminus:offsetQsminus+T-1) - x(offsetQsplus:offsetQsplus+T-1);
            sol.Inj(:,s) = [sol.Pb(:,s); sol.Pfc(:,s)-sol.Pel(:,s); sol.Qb(:,s); sol.Qfc(:,s)];
            % Get the voltage and current magnitudes
            Av = AVS{s}; Ai = AIS{s};
            sol.V(:,s) = Vbar(:,s) + Av*(sol.Inj(:,s)-Injbar(:,s)); sol.I(:,s) = Ibar(:,s) + Ai*(sol.Inj(:,s)-Injbar(:,s));
            
            % Check Loadflow
    
            for t=1:T
                Sinj = zeros(N,1);
                Sinj(1) = sol.ps(t,s) + 1i * sol.qs(t,s);
                Sinj(loads.indices(1)) = -10e3*loads.data(t,s)*1/Sb;
                Sinj(loads.indices(2)) = -5e3*loads.data(t,s)*1/Sb;
                Sinj(uncontrollable.indices(1)) = 10e3*uncontrollable.data{pv_num}(t,s)*1/Sb;
                Sinj(uncontrollable.indices(2)) = 10e3*uncontrollable.data{pv_num}(t,s)*1/Sb;
                Sinj(idxbat) = sol.Pb(t,s) + 1i * sol.Qb(t,s);
                Sinj(idxh) = sol.Pfc(t,s) - sol.Pel(t,s) + 1i * sol.Qfc(t,s);
                Vt = [1;sol.V(t:T:(N-2)*T+t,s)];
                [J,E,Snew,n_iter] = NR_polar_sol(Sinj,Vt,YY,E_0,idx,paramsNR);
                Errors(t,s) = norm(abs(E)-Vt);
                SLF(t,s) = real(Snew(1));
                Losses(t,s) = sum(real(Snew));
                ErrorsS(t,s) = abs(real(Snew(1))-sol.ps(t,s));
                Sbar(t:T:(N-1)*T+t,s) = Snew;
                pspbar(t,s) = max(0,real(Snew(1))); psnbar(t,s) = -min(0,real(Snew(1)));
                qspbar(t,s) = max(0,imag(Snew(1))); qsnbar(t,s) = -min(0, imag(Snew(1)));
                Vbar(t:T:(N-2)*T+t,s) = E(2:end); 
                Itemp = zeros(NL,1);
                for i=1:NL
                    idx1 = LP(i,1); idx2 = LP(i,2);
                    Itemp(i) = YYL(idx1,idx2)*(E(idx1)-E(idx2)) + YYT(idx1,idx2)*E(idx1);
                   % Itemp(i+NL) = YYL(idx2,idx1)*(E(idx2)-E(idx1)) +
                   % YYT(idx2,idx1)*E(idx2); Ideally also check second side
                end
                Ibar(t:T:(NL-1)*T+t,s) = abs(Itemp);
            end
        end
        Injbar = sol.Inj;
       disp('LF Errors:  '); disp(max(max(Errors))); disp(max(max(ErrorsS)));
       if max(max(ErrorsS)) <= tol 
            converged = 1;
        end
    end

% figure(1)
% hold on
% subplot(3,1,1)
% plot(erV(1:niter))
% subplot(3,1,2)
% plot(erS(1:niter))
% subplot(3,1,3)
% plot(erI(1:niter))
%% Displaying errors
disp('LF:  '); disp(max(Errors));
disp('Nb of subproblem inner iterations:'); disp(niter);
sol.ps = pspbar - psnbar; sol.qs = qspbar - qsnbar; 
sol.V = abs(Vbar); sol.I = abs(Ibar); 
sol.pb = sol.Pb;


%% Plotting

plotting = 0;
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
    phs = [sol.ph(1:end-1,scenario)];
    plot(time,phs*pb);
    title('Hydrogen Storage')
    xlabel('Time (hours)')
    ylabel('Pressure [bar]')
    
    
%% plot dispatch plan
figure(3)
hold on
for i = 1:9
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


%% Plot for all scenarios

for scenario=1:9
    figure(3+scenario)
    ActiveInj = cell(3,1); ActiveInj{1} = 'Battery'; ActiveInj{2} = 'Fuel Cell'; ActiveInj{3} =  'Electrolyzer';
    subplot(3,3,1)
    pbs = sol.Pb(:,scenario);
    psum = pbs*Sb/1000;
    plot(time,pbs*Sb/1000)
    title(ActiveInj(1))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    subplot(3,3,2)
    hold on
    pfcs = sol.Pfc(:,scenario);
    psum=psum + pfcs*Sb/1000;
    plot(time,pfcs*Sb/1000)
    title(ActiveInj(2))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    subplot(3,3,3)
    hold on
    pels = sol.Pel(:,scenario);
    psum=psum - pels*Sb/1000;
    plot(time,pels*Sb/1000)
    title(ActiveInj(3))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    
    for i=1:length(loads.indices)
        subplot(3,3,3+i)
        if i==1
            f=10e3;
        else
            f=5e3;
        end
        L = f*loads.data;
        Lds = L(1:T,scenario);
        psum=psum - Lds/1000;
        plot(time,-Lds/1000)
        xlabel('Time (hours)')
        title('Load')
        ylabel('Power [kW]')
    end
    for i=1:length(uncontrollable.indices)
        subplot(3,3,5+i)
        %U = uncontrollable.data{i};
        U = 10e3*uncontrollable.data{pv_num};
        Us = U(1:T);
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
    plot(time, value(dispatch)*Sb/1000)
    hold off
    title('Slack Power')
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    
    subplot(3,3,9)
    Ebs = [sol.Eb(1:end-1,scenario)];
    plot(time,Ebs*Sb/1000)
    title('Battery Energy')
    xlabel('Time (hours)')
    ylabel('Energy [kWh]')

end 
end  

if plotting == 2

    figure(2)
    for scenario=1:S
        subplot(3,3,scenario)
        hold on
        pbs = sol.Pb(:,scenario);
        psum = pbs*Sb/1000;
        plot(time,pbs*Sb/1000)
        title('Battery')
        xlabel('Time (hours)')
        ylabel('Power [kW]')
        hold off
    end 

    figure(3)
    for scenario=1:S
        subplot(3,3,scenario)
        hold on
        pfcs = sol.Pfc(:,scenario);
        psum=psum + pfcs*Sb/1000;
        plot(time,pfcs*Sb/1000)
        title('Fuel cell')
        xlabel('Time (hours)')
        ylabel('Power [kW]')
        hold off
    end 
    
    figure(4)
    for scenario=1:S
        subplot(3,3,scenario)
        hold on
        pels = sol.Pel(:,scenario);
        psum=psum - pels*Sb/1000;
        plot(time,pels*Sb/1000)
        title('Electrolyser')
        xlabel('Time (hours)')
        ylabel('Power [kW]')
        hold off
    end 
    
    figure(5)
    for scenario=1:S
        subplot(3,3,scenario)
        hold on
        Ebs = [sol.Eb(1:end-1,scenario)];
        plot(time,Ebs*Sb/1000)
        title('Battery Energy')
        xlabel('Time (hours)')
        ylabel('Energy [kWh]')
        hold off
    end 

figure(6)
for i = 1:S
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
    hold off
end

figure(7)
    for scenario=1:S
        subplot(3,3,scenario)
        hold on
        phs = [sol.Eh(1:end-1,scenario)];
        plot(time,sol.Eh(1:end-1,scenario)*Sb/1000);
        title('Hydrogen Storage')
        xlabel('Time (hours)')
        ylabel('Energy [kWh]')
        hold off
    end 

scenario=1;
figure(8)
    ActiveInj = cell(3,1); ActiveInj{1} = 'Battery'; ActiveInj{2} = 'Fuel Cell'; ActiveInj{3} =  'Electrolyzer';
    subplot(3,3,1)
    pbs = sol.Pb(:,scenario);
    psum = pbs*Sb/1000;
    plot(time,pbs*Sb/1000)
    title(ActiveInj(1))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    subplot(3,3,2)
    hold on
    pfcs = sol.Pfc(:,scenario);
    psum=psum + pfcs*Sb/1000;
    plot(time,pfcs*Sb/1000)
    title(ActiveInj(2))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    subplot(3,3,3)
    hold on
    pels = sol.Pel(:,scenario);
    psum=psum - pels*Sb/1000;
    plot(time,pels*Sb/1000)
    title(ActiveInj(3))
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    
    for i=1:length(loads.indices)
        subplot(3,3,3+i)
        %L = loads.data{i};
        if i==1
            f=10e3;
        else
            f=5e3;
        end
        L = f*loads.data;
        Lds = L(1:T,scenario);
        psum=psum - Lds/1000;
        plot(time,-Lds/1000)
        xlabel('Time (hours)')
        title('Load')
        ylabel('Power [kW]')
    end
    for i=1:length(uncontrollable.indices)
        subplot(3,3,5+i)
        %U = uncontrollable.data{i};
        U = 10e3*uncontrollable.data{pv_num};
        Us = U(1:T);
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
    plot(time, value(dispatch)*Sb/1000)
    hold off
    title('Slack Power')
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    
    subplot(3,3,9)
    Ebs = [sol.Eb(1:end-1,scenario)];
    plot(time,Ebs*Sb/1000)
    title('Battery Energy')
    xlabel('Time (hours)')
    ylabel('Energy [kWh]')

    
end

if plotting == 3
    figure
    sum_without_controllables = 20e3*uncontrollable.data{1}(:, 1:S) - 10e3*loads.data(:, 1:S) - 5e3*loads.data(:, 1:S);
    plot(time,sum_without_controllables)
    title('Uncontrollables')
    xlabel('Time (hours)')
    ylabel('Energy [Wh]')

    figure
    plot(time,sol.Pb*Sb/1000)
    title('Battery')
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off


    figure
    plot(time,sol.Pfc*Sb/1000)
    title('Fuel cell')
    xlabel('Time (hours)')
    ylabel('Power [kW]')
    hold off
    
    
    figure
    pels = sol.Pel;
    plot(time,pels*Sb/1000)
    title('Electrolyser')
    xlabel('Time (hours)')
    ylabel('Power [kW]')
        
    
    figure
    plot(time,sol.Eb(1:end-1,:)*Sb/1000)
    title('Battery Energy')
    xlabel('Time (hours)')
    ylabel('Energy [kWh]')
        

    figure
    plot(time, value(dispatch)*Sb/1000)
    hold on
    plot(time, sol.ps*Sb/1000)
    hold off
    title('Slack Power')
    xlabel('Time (hours)')
    ylabel('Power [kW]')

    legend('dispatch plan', 'Scenario 1', 'Scenario 2', 'Scenario 3', 'Scenario 4')
    hold off
    
    
    figure
    plot(time,sol.Eh(1:end-1,:)*Sb/1000);
    title('Hydrogen Storage')
    xlabel('Time (hours)')
    ylabel('Energy [kWh]')
        
    

end


fprintf('Ebmax: %.6f\n', sol.Ebmax*Sb);
fprintf('Sbmax: %.6f\n', sol.Sbmax*Sb);
fprintf('Ehmax: %.6f\n', sol.Ehmax*Sb);

fprintf('PFCmax: %.6f\n', sol.PFCmax*Sb);
fprintf('PELmax: %.6f\n', sol.PELmax*Sb);

end