%%% Main file: Call the Simplified Subproblem and solve for one scenario %%
close all;
load MicroGrid.mat
load Resourceparameters.mat

load microloadsreal.mat
loads.data = Loads_Scenarios; loads.indices = load_indices;
uncontrollable.data = UncGen_Scenarios; uncontrollable.indices = uncon_indices;
state.SOCH = 0.5; state.SOCB = 0.5; % Initial battery level, initial hydrogen pressure
scenario = 1;
T = 24*4; S = 10; Sb = 600000; Vb = 400;

% Construct aggregated load as dispatch
avload = zeros(T,1); avunc = zeros(T,1);
minload = zeros(T,1); minunc = zeros(T,1);
maxload = zeros(T,1); maxunc = zeros(T,1);
aggload = zeros(T,S);
for i=1:length(loads.indices)
    avload = avload + mean(loads.data{i}(1:T,1:S),2);
    mini = zeros(T,1); maxi = zeros(T,1);
    aggload = aggload + loads.data{i}(1:T,1:S);
    for t= 1:T
        mini(t) = min(loads.data{i}(t,1:S));
        maxi(t) = max(loads.data{i}(t,1:S));
    end
    minload = minload + mini;
    maxload = maxload + maxi;
end
for i=1:length(uncontrollable.indices)
    avunc = avunc + mean(uncontrollable.data{i}(1:T,1:S),2);
    mini = zeros(T,1); maxi = zeros(T,1);
    aggload = aggload - uncontrollable.data{i}(1:T,1:S);
    for t= 1:T
        mini(t) = min(uncontrollable.data{i}(t,1:S));
        maxi(t) = max(uncontrollable.data{i}(t,1:S));
    end
    minunc = minunc + mini;
    maxunc = maxunc + maxi;
end
aggLoad = 1/Sb * (-avunc + avload);
LoadTop = 1/Sb * (-minunc + maxload);
LoadBot = 1/Sb * (-maxunc + minload);


[sol, iterationcost] = SimplifiedSubproblem(loads,uncontrollable,scenario, grid, resources, costs, state, T, aggLoad);