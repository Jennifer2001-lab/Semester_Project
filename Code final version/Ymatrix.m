function [YY, YYL, YYT, Ib, Ampacities] = Ymatrix(Text,Ab,Vb)

%% Admittance matrix computation

% disp('                                                 ')
% disp('*************************************************')
% disp('*This script computes the [Y] admittance matrix.*')
% disp('*************************************************')
% disp('                                                 ')

%----------------------------------------------------------------------------------------------%
%% First step: Line parameters
%%%The information of the line parameters are loaded using a txt file 
%%%The lines of the txt file should have the following format:
%%%Bus Bus  R(Ohm/Km)   X(Ohm/Km)   B(S/Km) length(Km)
% disp('STEP 1: Loading line parameters....')
% disp('                                   ')

%linedata = load ('linedata_microgrid.txt'); %%% here add the name of the txt file with the line parameters
linedata = load (Text);

n_lines=length(linedata(:,1));                          %%%number of lines
n_nodes=max(max(linedata(:,1)),max(linedata(:,2))); %%%number of nodes

% disp(['The network consists of ',num2str(n_nodes), ' nodes and ',num2str(n_lines), ' lines.'])
% disp('                                                                                      ')

%----------------------------------------------------------------------------------------------%
%% Second step: Base values of the network
%%%The base values for the power, the voltage need to be defined here in
%%%order to compute the [Y] matrix in p.u.
%%%The base values for the current and the admmitance are directly computed

% disp('STEP 2: Transforming line parameters in p.u.....')
% disp('                                                ')

%Ab=100e3;    %%%base value for the power in VA
%Ab=1e5;    %%%base value for the power in VA
%Vb=400;  %%%base value for the line-to-line voltage in V

Ib=Ab/(Vb.*sqrt(3));
Yb=Ab/Vb^2;
% Zb=3.46112;
% Yb=1/Zb;

time=tic;

%%%The line parameters are expressed here in p.u.
line_lengths=linedata(:,6);
linedata(:,3)=(linedata(:,3).*Yb).*line_lengths;
linedata(:,4)=(linedata(:,4).*Yb).*line_lengths;
linedata(:,5)=(linedata(:,5)./Yb).*line_lengths;

%----------------------------------------------------------------------------------------------%
%% Third step: Building the primitive branch and shunt admittance matrices

% disp('STEP 3: Building the network primitive branch and shunt admittance matrices....')
% disp('                                                            ')

y_i_ih=zeros(n_nodes,n_nodes);
y_i=zeros(1,n_nodes); 

y_ih=1./(linedata(:,3)+1i*linedata(:,4)); %%%the longitudinal admittance of each line
absy_ih=abs(y_ih);
angley_ih=angle(y_ih);
g_ih=real(y_ih);
b_ih=imag(y_ih);


for k=1:n_lines
    y_i_ih(linedata(k,1),linedata(k,2))=linedata(k,5)*1i/2;
    y_i_ih(linedata(k,2),linedata(k,1))=linedata(k,5)*1i/2;
end

for k=1:n_nodes
    y_i(k)=sum(y_i_ih(k,:)); %%%the transversal admittance elements
end

g_i_ih=real(y_i_ih);
b_i_ih=imag(y_i_ih);

% Extra
for k=1:n_lines
    % This is not 100% Correct I am assuming that both shunt elements are
    % equal (reasonable...)
    y_i_cur(k)=y_i(linedata(k,1)) ;%+ y_i(linedata(k,2)); 
end

YL = diag(y_ih);
YT = diag(y_i);
YT_cur = diag(y_i_cur);

%----------------------------------------------------------------------------------------------%
%% Fourth step: compute incidence matrix of the network

% disp('STEP 4: Building the network incidence matrix....')
% disp('                                                 ')


%%% The branch-to-node incidence matrix A is computed here
%%% A is of size (number of branches, number of nodes)
%%% A_ij=0 if branch i is not connected to node j
%%% A_ij=1 if current in branch i is directed away from node j
%%% A_ij=-1 if current in branch i is directed towards node j

A=zeros(n_lines,n_nodes); 

for k=1:n_lines
    A(k,linedata(k,1))=1;
    A(k,linedata(k,2))=-1;
end


%% Fifth step: compute [Y] matrix and visualize network connectivity

% disp('STEP 5: Computing the [Y] matrix....')
% disp('                                                 ')

Y=A.'*YL*A + YT;

% For Current Sensitivites
YY = Y;
YYL = -(YY - diag(diag(YY)));
YYT = -(A.'*YT_cur*A - diag(diag(A.'*YT_cur*A)));
Ampacities = (-(A.'*diag(linedata(:,7))*A - diag(diag(A.'*diag(linedata(:,7))*A))))/Ib;


% %%%visualization of the network graph as described by Y matrix
% Y2=Y;
% Y2(logical(eye(size(Y,1))))=0;
% gObj = biograph(abs(Y2));
% view(gObj);
% 
% figure;
% spy(Y)
% axis([1 size(Y,1) 1 size(Y,1)])
% title('Sparsity Pattern of [Y]')

end
