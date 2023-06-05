function [J,E,S,n_iter] = NR_polar(S_star,E_star,Y,E_0,idx,Parameters)
%
% INPUT
% - Y           nodal admittance matrix
% - S_star      given complex powers (active/reactive powers)
% - E_star      given voltage magnitudes
% - E_0         initial voltages (phasors)
% - idx.slack   index of the slack bus
% - idx.pq      indices of the PQ buses
% - idx.pv      indices of the PV buses
% - Parameters.tol         tolerance for convergence criterion
% - Parameters.n_max       maximum number of iterations
%
% OUTPUT
% - E           solution voltages (phasors)
% - J           Jacobian at the solution
% - n_iter      number of iterations

Yabs = abs(Y);
Yarg = angle(Y);
n_nodes = length(E_0);

% Initialization
Eabs = abs(E_0);
Eabs(idx.pv) = E_star(idx.pv); % Set given values for PV nodes
Earg = angle(E_0);
J = [];

for k=1:Parameters.n_max
    n_iter = k;
    
    % Compute nodal voltages/currents/power
    E = complex(Eabs.*cos(Earg),Eabs.*sin(Earg));
    I = Y*E;
    S = E.*conj(I);
    
    %% Mismatch calculation
    
    % Compute the mismatches for the entire network.
    dS = S_star-S;
    dP = real(dS);
    dQ = imag(dS);
    
    % Keep only the relevant mismatches.
    dP(idx.slack) = [];
    dQ(sort([idx.pv;idx.slack])) = [];
    
    dF = [dP;dQ]; % mismatch of the power flow equations
    
    %% Convergence check
    
    if(max(abs(dF))<Parameters.tol)
        % disp('NR algorithm has converged to a solution!');
        break;
    elseif(k==Parameters.n_max)
        disp('NR algorithm reached the maximum number of iterations!');
    end
    
    %% Jacobian construction
    
    % For the sake of simplicity, the blocks of J are constructed
    % for the whole network (i.e., with size n_nodes x n_nodes).
    % The unnecessary rows/columns are removed subsequently
    
    % Extract magnitude/angle
    Eabs = abs(E);
    Earg = angle(E);
    
    % Initialization
    J_PE = zeros(n_nodes,n_nodes); % derivative: P versus E_abs
    J_PT = zeros(n_nodes,n_nodes); % derivative: P versus E_arg (theta)
    J_QE = zeros(n_nodes,n_nodes); % derivative: Q versus E_abs
    J_QT = zeros(n_nodes,n_nodes); % derivative: Q versus E_arg (theta)
    
    % Construction
    for i=1:n_nodes
        
        % Diagonal elements (terms outside the sum)
        J_PE(i,i) =  2*Yabs(i,i)*Eabs(i)*cos(Yarg(i,i));
        J_QE(i,i) = -2*Yabs(i,i)*Eabs(i)*sin(Yarg(i,i));
        
        for j=1:n_nodes
            if i ~= j
                % Diagonal elements (terms inside the sum)
                J_PE(i,i) = J_PE(i,i) + Yabs(i,j)*Eabs(j)*cos(Earg(i)-Earg(j)-Yarg(i,j));
                J_QE(i,i) = J_QE(i,i) + Yabs(i,j)*Eabs(j)*sin(Earg(i)-Earg(j)-Yarg(i,j));
                J_PT(i,i) = J_PT(i,i) - Eabs(i)*Yabs(i,j)*Eabs(j)*sin(Earg(i)-Earg(j)-Yarg(i,j));
                J_QT(i,i) = J_QT(i,i) + Eabs(i)*Yabs(i,j)*Eabs(j)*cos(Earg(i)-Earg(j)-Yarg(i,j));

                % Offdiagonal elements
%                 J_PE(i,j) =  Y_abs(i,j)*E_abs(i)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));
%                 J_QE(i,j) =  Y_abs(i,j)*E_abs(i)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
%                 J_PT(i,j) =  Y_abs(i,j)*E_abs(i)*E_abs(j)*sin(E_arg(i)-E_arg(j)-Y_arg(i,j));
%                 J_QT(i,j) = -Y_abs(i,j)*E_abs(i)*E_abs(j)*cos(E_arg(i)-E_arg(j)-Y_arg(i,j));
                
                J_PE(i,j) = J_PE(i,j)+Yabs(i,j)*Eabs(i)*cos(Earg(i)-Earg(j)-Yarg(i,j)); 
                J_QE(i,j) = J_QE(i,j)+Yabs(i,j)*Eabs(i)*sin(Earg(i)-Earg(j)-Yarg(i,j)); 
                J_PT(i,j) = J_PT(i,j)+ Yabs(i,j)*Eabs(i)*Eabs(j)*sin(Earg(i)-Earg(j)-Yarg(i,j));
                J_QT(i,j) = J_QT(i,j)-Yabs(i,j)*Eabs(i)*Eabs(j)*cos(Earg(i)-Earg(j)-Yarg(i,j)); 

            end
        end
    end
    
    % Remove extra rows (i.e., unnecessary equations)
    % slack bus: P & Q, PV buses: Q
    
    J_PE(idx.slack,:) = [];
    J_PT(idx.slack,:) = [];
    
    J_QE([idx.pv;idx.slack],:) = [];
    J_QT([idx.pv;idx.slack],:) = [];
    
    % Remove extra columns (i.e., variables)
    % slack bus: E_abs & E_arg, PV nodes: E_abs
    
    J_PE(:,[idx.pv;idx.slack]) = [];
    J_QE(:,[idx.pv;idx.slack]) = [];
    
    J_PT(:,idx.slack) = [];
    J_QT(:,idx.slack) = [];
    
    % Combination
    J = [J_PE,J_PT;J_QE,J_QT];
    
    %% Solution update
    
    % Solve
    dx = J \ dF;
    
    % Reconstruct the solution
    
    dEabs = zeros(length(Eabs),1);
    dEabs(idx.pq,1) = dx(1:length(idx.pq));
    
    dEarg = zeros(length(Earg),1);
    dEarg(sort([idx.pq;idx.pv]),1) = dx((length(idx.pq)+1):end);
    
    % Update
    Eabs = Eabs + dEabs;
    Earg = Earg + dEarg;
end

E = Eabs .* exp(1i*Earg);

end