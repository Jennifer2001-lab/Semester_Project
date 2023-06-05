function createResourceParameters()

load('MicroGrid.mat', 'grid' )
bv= grid.basevalues;
Vb=bv(1); Sb=bv(2); pb=bv(5); Ib = Sb/(sqrt(3)*Vb);

%%% Parameters of the resources
% Physical Parameters
basic.Mh2 = 2*1.00794; % g/mol
basic.F = 9.65*10^4; % Coulomb per mol
basic.Rgas = 4124; % Specific gas constant (J/kgK) from FC1

% Fuel Cell Parameters
resources.PFCmax =18000/Sb; resources.PFCmin = 0.1*resources.PFCmax; % Max/Min Operating point of the fuel cell 
%c etait 0 avant
resources.QFCmax = resources.PFCmax * 0.5; resources.QFCmin = -resources.QFCmax;
resources.kh_fc = 1/0.75*1/33.6; % simple model for hydrogem consumption

% Electrolyzer Parameters
resources.PELmax = 4500/Sb; resources.PELmin = 0.1*resources.PELmax; % Electrolyzer power limits 
%c etait 0 avant
resources.kh_el = 0.75*1/33.6; % simple model for hydrogen production

% Hydrogen Storage Parameters
Vtank = 1; % H2 tank volume in cubic meters
Tt = 293; % Assumed tank temperature
resources.phmax = 1*30/pb; resources.phmin = 2/pb; % Max hydrogen pressure in bar
resources.Ktank = (10^-5) * basic.Rgas * Tt/Vtank *Ib/pb; % Hydrogen production to pressure (in per unit)

% Battery Parameters
f=1; % a checker, plus de storgae mieux pour dipatch
resources.Ebmax = f*25000/Sb; resources.Ebmin = resources.Ebmax/5; %Wh
resources.Pbmax = resources.Ebmax; resources.Pbmin = -resources.Ebmax; %W
resources.Qbmax = resources.Pbmax*0.5; resources.Qbmin = resources.Pbmin * 0.5;

% Grid Parameters
resources.Pgmax = 35000/Sb; resources.Pgmin = -35000/Sb; % Maximum grid power
resources.theta_max = pi/6;

% Cost Parameters price in CHF/kWh
costs.ccurt = 1e6; costs.cshed = 1e5;
costs.closs = Sb/1000*1;
costs.cslack = 10*Sb/10;
costs.cfc = 5*0.05*Sb/1000; 
costs.cel = 5*0.05*Sb/1000; % price in CHF/kWh
costs.csfc = 1;% * Sb/1000; 
costs.csel = 1;% * Sb/1000;
costs.cbat = 2 * 0.05*Sb/1000; % Price of cycling energy through the battery CHF/kWh
costs.el_plus = 1;
costs.el_min= 1;
% NR Parameters
paramsNR.n_max = 1000; paramsNR.tol = 1e-6; 

save Resourceparameters.mat resources costs paramsNR basic

end