%%% Test solving subproblem %%%
close all;
clear;

load MicroGrid.mat
load Resourceparameters.mat 
load microloadsreal.mat load_indices uncon_indices
load Season_scenarios_weeks_Jennifer.mat
load Load_scenarios_weeks_Jennifer.mat 

% addpath(genpath('C:\Users\admin\Documents\rgupta\Yalmip\'))
% addpath('C:\Program Files\Mosek\10.0\toolbox\r2017a')
% addpath('C:\gurobi1001\win64\matlab\')

%% Load the data
loads.data = load_scenarios;
loads.indices = load_indices;
uncontrollable.data = pv_scenarios;
uncontrollable.indices = uncon_indices;

el_price = elpricenextday; 

T=24*7; % number of timesteps
S=8; % number of scenarios, max 8

pv_num = 1; subprob_number = 1;
% create an artificial master solution
Sb = 600000;
% master_sol.Ebmax = 25000/Sb; master_sol.Sbmax = master_sol.Ebmax;
% master_sol.Ehmax = 200000/Sb; master_sol.PFCmax = 20000/Sb; master_sol.PELmax = 10000/Sb;
% master_sol.fc_bin(1,1) = 1; master_sol.el_bin(1,1) = 1;

master_sol.Ebmax = 0; master_sol.Sbmax = 0;
master_sol.Ehmax = 0; master_sol.PFCmax = 0; master_sol.PELmax = 0;
master_sol.fc_bin(1,1) = 0; master_sol.el_bin(1,1) = 0;
%% Solve a subproblem

[subproblem_objective, subprob_solution_temp, dual_solution_sub] = solve_subproblem(loads, uncontrollable, grid, resources, costs, el_price, master_sol, T, S, pv_num, subprob_number)
