function createGridData(lineparamstext,Sb,Vb,lpsize)
%createGridData('Copy_of_linedata_microgrid_new.txt',600000,400,[7,inf])
Zb = Vb^2/Sb; Yb = 1/Zb; pb = Sb/1000;

formatSpec = '%f'; fileID = fopen(lineparamstext,'r');
LP = fscanf(fileID, formatSpec, lpsize)';

NL = size(LP,1); N = max(max(LP(:,1:2))); 
NC = 5; % Number of controllable generators considered (active and reactive)

idx_slack=1; idx_pv=[]; idx_pq=[2:N];
idx.slack=idx_slack; idx.pv=idx_pv; idx.pq=idx_pq;

[YY, YYL, YYT, Ib, Ampacities] = Ymatrix(lineparamstext,Sb,Vb);

% Create R,X,B vectors for AR OPF
R=zeros(NL,1); X=zeros(NL,1); B=zeros(NL,1);
for i=1:NL
   idx1=LP(i,1); idx2=LP(i,2);
   Z = 1./YYL(idx1,idx2);
   R(i)= real(Z); X(i) = imag(Z);
   B(i)=imag(YYT(idx1,idx2));
end

% Create connection matrices and node to line linking matrices
A = zeros(NL); Ab = zeros(NL,N); Ae = zeros(NL,N);
for i=1:NL
    branchstart = LP(i,1); branchend = LP(i,2);
    Ab(i,branchstart) = 1; Ae(i,branchend) = 1;
    starting = find(LP(:,1)==branchend);
    for j=1:length(starting)
        startbranch = starting(j);
        A(i,startbranch) = 1;
    end
end
Abe = Ae-Ab;

% Save all parameters in a structure
grid.lineparameters = LP;
grid.basevalues = [Vb, Sb, Zb, Yb, pb];
grid.lines = NL; grid.nodes = N; grid.controllable = NC;
grid.node_indices = idx;
grid.Admittance = cell(5,1);
grid.Admittance{1} = YY; grid.Admittance{2} = YYL; grid.Admittance{3} = YYT;
grid.Admittance{4} = Ib; grid.Admittance{5} = Ampacities;
grid.connections = cell(4,1);
grid.connections{1} =A; grid.connections{2}=Ab; grid.connections{3}=Ae; grid.connections{4}=Abe;
grid.resistances = [R, X, B];

%save MicroGrid.mat grid
save MicroGrid_new.mat grid

end