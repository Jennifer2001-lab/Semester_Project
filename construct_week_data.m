clear 

load Season_scenarios_Jennifer.mat load_scenarios
load Season_scenarios_weeks_Jennifer.mat
load SimpleScenarios.mat elpricenextday

% PV data
uncontrollable.data = pv_scenarios;
% {1,1}, {1,2}, {1,3}, {1,4}
size(pv_scenarios{1});
% 8 scenarios for each
% hour resolution

% Load data
loads.data = load_scenarios;
size(load_scenarios{1});
size(load_scenarios{2});

load_weekday_res_15min = load_scenarios{1};
load_weekend_res_15min = load_scenarios{2};

load_weekday = zeros(24, 20);
load_weekend = zeros(24, 20);

% Change resolution
for i=0:(96/4 - 1)
    load_weekday(1+i, :) = load_weekday_res_15min( 1 + 4*i, :);
    load_weekend(1+i, :) = load_weekend_res_15min(1 + 4*i, :);
end

size(load_weekday);
size(load_weekend);

% Concatenate data
load_scenarios_week = [load_weekday; load_weekday; load_weekday; load_weekday; load_weekday; load_weekend; load_weekend];
size(load_scenarios_week);

load_scenarios = load_scenarios_week;


size(elpricenextday);
elprice_res_1hour = zeros(1, 24);
% Change resolution
for i=0:(96/4 - 1)
    elprice_res_1hour(1+i) = elpricenextday(1 + 4*i);
end
size(elprice_res_1hour);

% Concatenate data
elpricenextday = repmat(elprice_res_1hour, 1, 7);
size(elpricenextday);

save Load_scenarios_weeks_Jennifer.mat load_scenarios elpricenextday

