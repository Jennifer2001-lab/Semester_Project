function [K_p,K_com_p,K_q,K_comp_q]=Coeffs_Voltage_Alpha(Y,S0,E,Res_nodes_no_slack,slack,nph,alphap,alphaq,E0)

% This compute the exact voltage sensitivity coefficients
% Inputs
% YY                    admittance matrix
% S0                    nodal apparent power injections at nominal voltage
% E                     nodal voltages (state of the grid)
% Res_nodes_no_slack    a vector with "monophase" indexing without with the
% location of the ressources to which the user wants to compute the voltage
% sensitivity coefficients. Usually the slack node is the 1st node
% therefore if orginially node 6 has a resources then here one should input
% 5 in the vector. The output will be a nphN x nph matrix where N is the number of nodes (monophase)
% and nph are the number of phases, contataning the partial derivatives.
% slack                 monophase index of slack
% nph                   number of phases
% alphap and alphaq     vectors contataining injection voltage dependency
% E0                    constant for voltage dependency
% Outputs
% K_p                   Magnitude voltage sensitivity coefficients w.r.t P
% K_com_p               Complex voltage sensitivity coefficients w.r.t P
% K_q                   Magnitude voltage sensitivity coefficients w.r.t Q
% K_com_q               Complex voltage sensitivity coefficients w.r.t Q
% Sizes of all outputs are nphN x (number_of_ressources x nph)

a=(slack-1)*nph+1;
b=(slack-1)*nph+nph;

L=size(Y);
Enew=E(1:L);
Enew(a:b)=[];
absEnew = abs(Enew);  % Added
Ynew=Y(1:L,1:L);
Ynew(a:b,:)=[];
Ynew(:,a:b)=[];
Snew=S0(1:L);
Snew(a:b)=[];
alphap = alphap(1:L);
alphap(a:b) = [];
alphaq = alphaq(1:L);
alphaq(a:b) = [];

Pnew = real(Snew); Qnew = imag(Snew);

%MATRICES A,B,C,D contain the multipliers of the real and imaginary partial
%derivatives
F=Y*E;
Fn=F(1:L);
Fn(a:b)=[];

G=[];
for i=1:1:size(Enew)
    G(i,:)=conj(Enew(i))*Ynew(i,:);
end

% Computing the new system with the equations already divided.
finA=[];
for i=1:1:size(Enew)
    for j=1:1:size(Enew)
        if i==j
            % Real Equations
            Const_Real_equ = alphap(i)*(absEnew(i)^(alphap(i)-2))*Pnew(i) - 1i*alphaq(i)*(absEnew(i)^(alphaq(i)-2))*Qnew(i); 
            finA(i,j)=(real(Fn(i))+real(G(i,j))) - real(Enew(i))*real(Const_Real_equ) ; % Add here new term
            finA(i,size(Enew,1)+j)=(imag(Fn(i))-imag(G(i,j))) - imag(Enew(i))*real(Const_Real_equ)  ; % Add here new term 
            
            % Imaginary Equations
            Const_Imag_equ = alphap(i)*(absEnew(i)^(alphap(i)-2))*Pnew(i) - 1i*alphaq(i)*(absEnew(i)^(alphaq(i)-2))*Qnew(i); 
            finA(size(Enew,1)+i,j)=(imag(Fn(i))+imag(G(i,j))) - real(Enew(i))*imag(Const_Imag_equ) ; % Add here new term
            finA(size(Enew,1)+i,size(Enew,1)+j)=(real(G(i,j))-real(Fn(i))) - imag(Enew(i))*imag(Const_Imag_equ) ; % Add here new term 
        else
            finA(i,j)=real(G(i,j));
            finA(size(Enew,1)+i,j)=imag(G(i,j));
            finA(i,size(Enew,1)+j)=-imag(G(i,j));
            finA(size(Enew,1)+i,size(Enew,1)+j)=real(G(i,j));
        end
    end
end

%computing the coefficients for the active powers P
kk=1; K_p=[]; K_com_p=[];
for k=1:length(Res_nodes_no_slack)
     for l = 1:nph
        % Construction of the constants
        RHS=zeros(2*size(Enew,1),1);
        current_idx = (Res_nodes_no_slack(k)-1)*nph; % + l; EDIT
        RHS(current_idx)=(absEnew(current_idx)/E0(k+1))^alphap(current_idx); 
        
        % Solve the linear system
        y=linsolve(finA,RHS);
        K_com_p(:,kk)=y;
        % Computing the coefficients of delta_|E|/delta_P
        K_p(:,kk)=(1./abs(Enew)).*real( complex(y(1:size(Enew)),-y(size(Enew)+1:2*size(Enew))).*Enew);

        % Index Update
        kk=kk+1;
     end
end

%COMPUTING THE COEFFICIENTS FOR Q
ll=1; K_q=[]; K_comp_q=[];
for k=1:length(Res_nodes_no_slack)
     for l = 1:nph
        % Construction of the constants
        RHS=zeros(2*size(Enew,1),1);
        current_idx = (Res_nodes_no_slack(k)-1)*nph; %  + l; EDIT
        RHS(current_idx+(size(Enew,1)))=-(absEnew(current_idx)/E0(k+1))^alphaq(current_idx); 
        
        
        % Solve the linear system
        y=linsolve(finA,RHS);
        K_comp_q(:,ll)=y;
        %computing the coefficients of delta_|E|/delta_Q
        K_q(:,ll)=(1./abs(Enew)).*real( complex(y(1:size(Enew)),-y(size(Enew)+1:2*size(Enew))).*Enew);
        %end
        ll=ll+1;
     end
end

% ReConstruct Complex Voltage Coeffs
K_com_p = K_com_p(1:0.5*size(K_com_p,1),:) ...
                  + 1i*K_com_p(0.5*size(K_com_p,1)+1:end, :);
K_comp_q = K_comp_q(1:0.5*size(K_comp_q,1),:) ...
                  + 1i*K_comp_q(0.5*size(K_comp_q,1)+1:end, :); 
              
% Add Slacks
K_com_p = [K_com_p(1:(nph*(slack-1) + 1 - 1),:); zeros(nph,size(K_com_p,2)) ; K_com_p(nph*(slack-1)+1:end,:)];
K_comp_q = [K_comp_q(1:(nph*(slack-1) + 1 - 1),:); zeros(nph,size(K_comp_q,2)) ; K_comp_q(nph*(slack-1)+1:end,:)];

K_p = [K_p(1:(nph*(slack-1) + 1 - 1),:); zeros(nph,size(K_p,2)) ; K_p(nph*(slack-1)+1:end,:)];
K_q = [K_q(1:(nph*(slack-1) + 1 - 1),:); zeros(nph,size(K_q,2)) ; K_q(nph*(slack-1)+1:end,:)];
end