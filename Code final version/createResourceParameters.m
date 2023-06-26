function createResourceParameters()

load('MicroGrid.mat', 'grid' )
bv= grid.basevalues;
Vb=bv(1); Sb=bv(2); pb=bv(5); Ib = Sb/(sqrt(3)*Vb);

%%% Parameters of the resources
% Physical Parameters
basic.Mh2 = 2*1.00794; % g/mol
basic.F = 9.65*10^4; % Coulomb per mol
basic.R = 8.314; %J/(kg K)
basic.Rgas = 8.314/(1e-3 * basic.Mh2); % Specific gas constant (J/molK) from FC1
basic.LHV_h2 = 33300; %Wh/kg

% Fuel Cell Parameters
resources.PFCmax =18000/Sb; resources.PFCmin = 0.1*resources.PFCmax; % Max/Min Operating point of the fuel cell 
resources.QFCmax = resources.PFCmax * 0.5; resources.QFCmin = -resources.QFCmax;
resources.IFCmax = resources.PFCmax/0.5; resources.IFCmin = 0.05*resources.IFCmax;
load FC_params.mat
resources.pfc = pfc;
resources.pfc.Ns = 18*resources.pfc.Ns;
resources.kh_fc = 3.6*1/(2*basic.F)*basic.Mh2*1.05*resources.pfc.Ns; % Fuel consumption (kg H2 per Ah, including 5% loss factor)
resources.Ih_fc = 1/0.95 * resources.pfc.Ns;
resources.aveffFC = 0.65;

resources.etafc = 0.65; resources.etael = 0.26;
% Electrolyzer Parameters
resources.PELmax = 4500/Sb; resources.PELmin = 0.3*resources.PELmax; % Electrolyzer power limits
resources.IELmax = resources.PELmax/0.2; resources.IELmin = 0.1*resources.IELmax;
load EL_params.mat
resources.pel = pel;
resources.pel.Ns = 5*resources.pel.Ns;
resources.kh_el = 3.6*1/(2*basic.F)*basic.Mh2*0.95*resources.pel.Ns; % Fuel production (kg H2 per Ah, including 5% loss factor)
resources.Ih_el = 0.95 * resources.pel.Ns;
resources.aveffEL = 0.65;

% Hydrogen Storage Parameters
Vtank = 1.6; % H2 tank volume in cubic meters
Tt = 293; % Assumed tank temperature
resources.phmax = 1*30/pb; resources.phmin = 2/pb; % Max hydrogen pressure in bar
resources.Ktank = (10^-5) * basic.Rgas * Tt/Vtank *Ib/pb;
resources.Kstorage = basic.LHV_h2 * 3.6 * basic.Mh2 / ( 2*basic.F) * Ib/Sb;
resources.Ehmax = resources.phmax * pb * 1e5* basic.Mh2 * 1e-3 * basic.LHV_h2 * Vtank / (basic.R * Tt)/Sb; % Pa*g/mol * Wh/kg * m^3 * 1/K * kg * K/J
resources.Ehmin = resources.Ehmax/15;


% Battery Parameters
f=1;
resources.Ebmax = f*25000/Sb; resources.Ebmin = resources.Ebmax/5; %Wh
resources.Pbmax = resources.Ebmax; resources.Pbmin = -resources.Ebmax; %W
resources.Qbmax = resources.Pbmax*0.5; resources.Qbmin = resources.Pbmin * 0.5;

% Second Battery Parameters
resources.Eb2max =0; resources.Eb2min = 0.2 * resources.Eb2max;
resources.Pb2max = resources.Eb2max; resources.Pb2min = -resources.Pb2max;
resources.Qb2max = resources.Pb2max*0.5; resources.Qb2min = resources.Pb2min * 0.5;
% Grid Parameters
resources.Pgmax = 35000/Sb; resources.Pgmin = -35000/Sb; % Maximum grid power
resources.theta_max = pi/6;

% Cost Parameters price in CHF/kWh
costs.ccurt = 100 * Sb/1e6; % default, replace by 95 % of electricity cost 
costs.cshed = 5 * Sb/1e3; % VOLL, take 5 eur/kWh
costs.cviol = Sb/1000;
costs.closs = costs.ccurt;
costs.cslack = 10*Sb;
costs.pwafactor = 0.1*Sb/1000;
% For fuel cell, consider purchase cost of 40k and EOL at 150mV degradation
fccost = 4000; 
costs.cfc = 0.02/150 * fccost * 1/resources.PFCmax; 
costs.crfc = 0.05/150 * fccost * 1/resources.PFCmax;
costs.csfc = 0.2/150 * fccost; % Quite optimistic based on Mardit -> multiply by 5?
% Similarly for electrolyzer (100k as price)
elcost = fccost; % should be 100000
costs.cel = 0.02/150 * elcost * 1/resources.PELmax; 
costs.crel = 0.05/150 * elcost * 1/resources.PELmax;
costs.csel = 0.2/150 * elcost; 
% For battery assume 2000 cycles @ 1C and purchase cost of 1000/kWh
costs.cbat = f * 1/resources.Pbmax * 25 * 200/2000; % Price of cycling energy through the battery CHF/kWh assuming 200CHF/kWh
costs.crbat = 1/2*costs.cbat;
costs.cbat2 = costs.cbat;

% NR Parameters
paramsNR.n_max = 1000; paramsNR.tol = 1e-6; 

%%% Compute (Optimal) Piecewise Affine Hyperplane Approximations for FC and EL
% Temporary fix fuel cell and electrolyzer temperature
T_fc = 333; T_el = T_fc;
Napprox = 200; epsilon =1;
%[start_indices,end_indices, points, coeffsA_fc, coeffsB_fc] = OptSCFit_Current(T_fc,Napprox,resources.pfc,epsilon);
%[start_indices,end_indices, points, coeffsA_el, coeffsB_el] = OptSCFit_Current_EL(T_el,Napprox, resources.pel,epsilon);
%resources.coeffsA_fc = coeffsA_fc * 1/(sqrt(3)*Vb); resources.coeffsB_fc = coeffsB_fc*1/Sb;
%resources.coeffsA_el = coeffsA_el * 1/(sqrt(3)*Vb); resources.coeffsB_el = coeffsB_el*1/Sb;

save Resourceparameters.mat resources costs paramsNR basic

end