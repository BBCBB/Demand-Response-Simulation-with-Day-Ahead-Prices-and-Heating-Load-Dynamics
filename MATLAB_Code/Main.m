clc
clear all
close all

%% --------- Paths ----------
inputFolder = fullfile(pwd, 'inputdata');   % where your CSVs live

%% --------- Load file name lists (user's helpers) ----------
sub_loadfilenames1
sub_loadfilenames2
sub_loadfilenames3

%% --------- Month view: day-ahead (DA) and temperature ----------
DA_record = [];
temp_table = readtable(fullfile(inputFolder, file_names2{1}));
temp_record = [];

% NOTE: original code loops 1:31 across file_names1; keep structure but safe
for i = 1:31
    DA_table = readtable(fullfile(inputFolder, file_names1{i}));
    x = size(DA_table,1);
    N_data_DA = x/15;                       % 15 zones
    DA_reshaped = reshape(DA_table.LBMP___MWHr_, 15, N_data_DA).';  % [hours x 15]
    timestamps_raw = reshape(datetime(DA_table.TimeStamp), 15, N_data_DA).';

    % append temperature (column 16) matched by nearest date
    DA_reshaped(:,16) = 0;
    for l = 1:N_data_DA
        [~, idx] = min(abs(datenum(timestamps_raw(l,1)) - datenum(datetime(temp_table.DATE))));
        DA_reshaped(l,16) = temp_table.HourlyDryBulbTemperature(idx);
    end
    DA_record = [DA_record; DA_reshaped];
end
% DA_record now: (#hours in month) x 16 (15 zones + temp)

% 3-hour peaks across the month (no changes except robust indexing)
if size(DA_record,1) >= 744
    threehour = zeros(248,15);
    k = 1;
    for ii = 1:3:744
        threehour(k,:) = max(DA_record(ii:ii+2,1:15), [], 1);
        k = k + 1;
    end
    time = 1:248;
    figure(1)
    plot(time, threehour(:,1:15))
    legend('CAPITL','CENTRL','DUNWOD','GENESE','H Q','HUD VL','LONGIL','MHK VL','MILLWD','N.Y.C.','NORTH','NPX','O H','PJM','WEST','Location','bestoutside')
    xlabel('3-hour periods over one month'), ylabel('Peak Price')
    sgtitle('3-hour peak loads over Jan 2019');

    figure(2)
    plot(time, threehour(:,10))
    legend('N.Y.C. over 31 days (3-hour intervals in order)')
    xlabel('3-hour periods over one month'), ylabel('Peak Price')
    sgtitle('3-hour peak loads over Jan 2019 - NYC');

    threehour1day = zeros(8,15);
    for i = 1:8
        threehour1day(i,:) = max(threehour(i:8:248,1:15), [], 1);
    end
    time = 1:8;
    figure(3)
    plot(time, threehour1day(:,1:15))
    legend('CAPITL','CENTRL','DUNWOD','GENESE','H Q','HUD VL','LONGIL','MHK VL','MILLWD','N.Y.C.','NORTH','NPX','O H','PJM','WEST','Location','bestoutside')
    xlabel('3-hour periods - 1-day peak'), ylabel('Peak Price')
    sgtitle('3-hour peak loads - day time peak');
end

%% --------- Extract target week (7 days, hourly) ----------
DA_record = [];
temp_table = readtable(fullfile(inputFolder, file_names2{1}));
temp_table = fillmissing(temp_table, 'linear');

for i = 1:7
    DA_table = readtable(fullfile(inputFolder, file_names3{i}));
    x = size(DA_table,1);
    N_data_DA = x/15;
    DA_reshaped = reshape(DA_table.LBMP___MWHr_, 15, N_data_DA).';  % [hours x 15]
    DA_reshaped(:,16) = 0;                                          % temperature col

    timestamps_raw = reshape(datetime(DA_table.TimeStamp), 15, N_data_DA).';
    for l = 1:N_data_DA
        [~, idx] = min(abs(datenum(timestamps_raw(l,1)) - datenum(datetime(temp_table.DATE))));
        DA_reshaped(l,16) = temp_table.HourlyDryBulbTemperature(idx);
    end
    DA_record = [DA_record; DA_reshaped];
end
% DA_record is 168 x 16 for the week

%% --------- Build per-minute series for one zone (N.Y.C. = col 10) + temperature ----------
% To keep shapes exact, repeat each hourly value 60 times (1 value per minute)
zone_col = 10;                      % N.Y.C.
hourly_price = DA_record(:, zone_col);   % 168x1
hourly_tempF = DA_record(:, 16);         % 168x1

% minute series (10080x1)
price_min = repelem(hourly_price, 60);
tempF_min  = repelem(hourly_tempF, 60);

% Fahrenheit to Celsius
tempC_min = (tempF_min - 32) * 5/9;

% Pack into "interpolatedData" as earlier code expects: [minutes x 2] = [price, tempC]
interpolatedData = [price_min, tempC_min];     % 10080 x 2

%% --------- Single-unit simulation: Basecase over the week (per-minute) ----------
% We simulate hour by hour; inside each hour we run 60 steps and store 60 samples.
% Keep system parameters as in your code.
Tr_record_basecase = [];   % store room temp per minute (for the week)
Q_record_basecase  = [];   % store load per minute (for the week)

Time = 60;                 % steps per hour (minutes)
Timestep = 1;              % "seconds" in your original notation; kept for consistency
N_timestep = 60;           % loop count per hour (do NOT go to +1)
Tset = 22;                 % deg C
TDB  = 1;                  % deg C

Prated = 5020; Area = 228; Vair = 228*5;
DensityAir = 1.225; Cp_air = 1005; Cv_air = 780;
C_funiture = 1000; V_furniture = 1000;
UA = 111; UAmass = 3924;

R1 = 1/UA;
R2 = 1/UAmass;
Ca = Vair*DensityAir*Cv_air;
Cm = V_furniture*C_funiture;

A = [-(1/(R2*Ca)+1/(R1*Ca))  1/(R2*Ca); 1/(R2*Cm)  -1/(R2*Cm)];
B = [1/(R1*Ca) 1/Ca ;0 0];

% Simulate per hour (168 hours)
x0 = [22; 22];
for hh = 1:168
    % Take the hour's outdoor temp as the last minute's value of that hour
    Tout_hourC = tempC_min((hh-1)*60 + (1:60));    % 60x1
    Q = zeros(1, Time);
    status = zeros(1, Time);

    % Initialize ON/OFF based on setpoint
    if     x0(1) > (Tset + TDB/2), status(1) = 0;
    elseif x0(1) < (Tset - TDB/2), status(1) = 1;
    else,  status(1) = 0;
    end
    Q0 = (status(1)==1) * (Prated*Timestep);  % W·s (kept from your code)
    u0 = [Tout_hourC(1); Q0];

    Troom_record = zeros(1, N_timestep);   % minute temps
    Tmass_record = zeros(1, N_timestep);

    for i = 2:N_timestep
        DT = A*x0 + B*u0;
        x1 = x0 + DT;

        % Hysteresis
        if     (x1(1) > (Tset + TDB/2)) && (status(i-1) == 1)
            status(i) = 0;
        elseif (x1(1) < (Tset - TDB/2)) && (status(i-1) == 0)
            status(i) = 1;
        else
            status(i) = status(i-1);
        end

        % Electrical input
        if status(i) == 1
            Q(i) = Prated * Timestep;
        else
            Q(i) = 0;
        end

        % advance
        x0 = x1;
        u0 = [Tout_hourC(i); Q(i)];
        Troom_record(i-1) = x1(1);
        Tmass_record(i-1) = x1(2);
    end

    % Append exactly 60 samples for this hour
    Tr_record_basecase = [Tr_record_basecase; Troom_record];  % 1x60 row appended → 168x60
    Q_record_basecase  = [Q_record_basecase;  Q];             % 1x60 row appended → 168x60
end

%% --------- Plots: Basecase (per day, per minute) ----------
color = ['k' 'y' 'c' 'b' 'r' 'm' 'g'];

for d = 1:7
    series = reshape(Q_record_basecase((d-1)*24+1:d*24, :).', [], 1);  % 1440x1
    n = numel(series);
    Time_hr = (0:n-1)/60;  % hours 0..24

    figure(20+d)
    set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
    plot(Time_hr, series, color(d))
    xlabel('Time (hour)'); ylabel('Load (arb. units)')
    sgtitle(['Load (per minute) on day ' num2str(d) ' of the week'])
end

% Weekly temperature, hourly (sample every 60th minute)
figure(30);
Tr_vec = reshape(Tr_record_basecase.', [], 1);  % 10080x1
Time_hr = 0:167;
set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
plot(Time_hr, Tr_vec(1:60:end), color(4))
xlabel('Time (hour)'); ylabel('Temperature (^oC)')
hold on
for i = 1:7
    xline(15 + (i-1)*24, 'r--', 'LineWidth', 2);
    xline(18 + (i-1)*24, 'r--', 'LineWidth', 2);
end
sgtitle('Hourly room temperature over the week (Basecase)')

%% --------- DR Case 1 (forced OFF during 3–6 pm) ----------
Tr_record_case1 = [];
Q_record_case1  = [];

x0 = [22;22];
for hh = 1:168
    day = floor((hh-1)/24) + 1;
    hour_of_day = mod(hh-1,24);  % 0..23
    force_off = (hour_of_day >= 15) && (hour_of_day < 18);  % 3–6 pm

    Tout_hourC = tempC_min((hh-1)*60 + (1:60));
    Q = zeros(1,Time); status = zeros(1,Time);

    if     x0(1) > (Tset + TDB/2), status(1)=0;
    elseif x0(1) < (Tset - TDB/2), status(1)=1;
    else,  status(1)=0;
    end
    if force_off, status(1)=0; end

    Q0 = (status(1)==1) * (Prated*Timestep);
    u0 = [Tout_hourC(1); Q0];
    Troom_record = zeros(1,N_timestep);

    for i=2:N_timestep
        DT = A*x0 + B*u0; x1 = x0 + DT;

        if force_off
            status(i) = 0;
        elseif (x1(1) > (Tset + TDB/2)) && (status(i-1)==1)
            status(i)=0;
        elseif (x1(1) < (Tset - TDB/2)) && (status(i-1)==0)
            status(i)=1;
        else
            status(i)=status(i-1);
        end

        if status(i)==1, Q(i)=Prated*Timestep; else, Q(i)=0; end
        x0=x1; u0=[Tout_hourC(i); Q(i)];
        Troom_record(i-1)=x1(1);
    end

    Tr_record_case1 = [Tr_record_case1; Troom_record];
    Q_record_case1  = [Q_record_case1;  Q];
end

% Case 1 daily plots (per minute)
for d = 1:7
    series = reshape(Q_record_case1((d-1)*24+1:d*24, :).', [], 1);
    n = numel(series); Time_hr = (0:n-1)/60;

    figure(42+d)
    set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
    plot(Time_hr, series, color(d))
    xlabel('Time (hour)'); ylabel('Load (arb. units)')
    sgtitle(['Per-minute load on day ' num2str(d) ' (Case 1)'])
end

% Case 1 weekly temperature, hourly
figure(50)
Tr_vec = reshape(Tr_record_case1.', [], 1);
Time_hr = 0:167;
set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
plot(Time_hr, Tr_vec(1:60:end), color(1))
xlabel('Time (hour)'); ylabel('Temperature (^oC)')
hold on
for i=1:7
    xline(15 + (i-1)*24, 'r--', 'LineWidth', 2);
    xline(18 + (i-1)*24, 'r--', 'LineWidth', 2);
end
sgtitle('Hourly room temperature over the week (Case 1)')

%% --------- DR Case 2 (lower setpoint 18 C during 3–6 pm) ----------
Tr_record_case2 = [];
Q_record_case2  = [];

x0 = [22;22];
for hh = 1:168
    hour_of_day = mod(hh-1,24);
    if (hour_of_day >= 15) && (hour_of_day < 18)
        Tset0 = 18;
    else
        Tset0 = Tset;
    end

    Tout_hourC = tempC_min((hh-1)*60 + (1:60));
    Q=zeros(1,Time); status=zeros(1,Time);

    if     x0(1) > (Tset0 + TDB/2), status(1)=0;
    elseif x0(1) < (Tset0 - TDB/2), status(1)=1;
    else,  status(1)=0;
    end
    Q0 = (status(1)==1)*(Prated*Timestep);
    u0 = [Tout_hourC(1); Q0];
    Troom_record = zeros(1,N_timestep);

    for i=2:N_timestep
        DT = A*x0 + B*u0; x1 = x0 + DT;

        if     (x1(1) > (Tset0 + TDB/2)) && (status(i-1)==1)
            status(i)=0;
        elseif (x1(1) < (Tset0 - TDB/2)) && (status(i-1)==0)
            status(i)=1;
        else
            status(i)=status(i-1);
        end
        if status(i)==1, Q(i)=Prated*Timestep; else, Q(i)=0; end
        x0=x1; u0=[Tout_hourC(i); Q(i)];
        Troom_record(i-1)=x1(1);
    end

    Tr_record_case2 = [Tr_record_case2; Troom_record];
    Q_record_case2  = [Q_record_case2;  Q];
end

% Case 2 daily plots (per minute)
for d = 1:7
    series = reshape(Q_record_case2((d-1)*24+1:d*24, :).', [], 1);
    n = numel(series); Time_hr = (0:n-1)/60;

    figure(64+d)
    set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
    plot(Time_hr, series, color(d))
    xlabel('Time (hour)'); ylabel('Load (arb. units)')
    sgtitle(['Per-minute load on day ' num2str(d) ' (Case 2)'])
end

% Case 2 weekly temperature, hourly
figure(70)
Tr_vec = reshape(Tr_record_case2.', [], 1);
Time_hr = 0:167;
set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
plot(Time_hr, Tr_vec(1:60:end), color(1))
xlabel('Time (hour)'); ylabel('Temperature (^oC)')
hold on
for i=1:7
    xline(15 + (i-1)*24, 'r--', 'LineWidth', 2);
    xline(18 + (i-1)*24, 'r--', 'LineWidth', 2);
end
sgtitle('Hourly room temperature over the week (Case 2)')

%% --------- Weekly cost (hourly) ----------
% All series are minute-based (10080). Convert to hourly: reshape 60x168 and average.
total_cost = zeros(1,3);

Qb = reshape(Q_record_basecase.', [], 1);    % 10080x1
Q1 = reshape(Q_record_case1.',   [], 1);
Q2 = reshape(Q_record_case2.',   [], 1);

load_per_hour_b = sum(reshape(Qb, 60, 168), 1)/60;  % 1x168
load_per_hour_1 = sum(reshape(Q1, 60, 168), 1)/60;
load_per_hour_2 = sum(reshape(Q2, 60, 168), 1)/60;

% price is hourly: DA_record(:,10)
cost_per_hour_b = (load_per_hour_b .* DA_record(:,10).')/1e6;  % scaling kept from your code
cost_per_hour_1 = (load_per_hour_1 .* DA_record(:,10).')/1e6;
cost_per_hour_2 = (load_per_hour_2 .* DA_record(:,10).')/1e6;

total_cost(1) = sum(cost_per_hour_b);
total_cost(2) = sum(cost_per_hour_1);
total_cost(3) = sum(cost_per_hour_2);

disp('Total weekly cost [base, case1, case2]:');
disp(total_cost);

%% --------- 50 Space Heating Units (per minute) ----------
MINUTES_WEEK = 7*24*60;     % 10080
ngbhd_load = zeros(MINUTES_WEEK, 3);

% randomize per your distributions
rng(42)
Prated50 = round(normrnd(5020, 0.03*5020, [50,1]), 0);
Area50   = round(normrnd(228,  0.04*228,  [50,1]), 0);
Vair50   = Area50.*round(normrnd(5,0.02*5,[50,1]),1);
Density50= round(normrnd(1.225,0.01*1.225,[50,1]),3);
Cp50     = round(normrnd(1005, 0.02*1005,[50,1]),0);
Cv50     = round(normrnd(780,  0.02*780, [50,1]),0);
Cf50     = round(normrnd(1000, 0.02*1000,[50,1]),0);
Vf50     = round(normrnd(1000, 0.02*1000,[50,1]),0);
UA50     = round(normrnd(111,  0.03*111, [50,1]),0);
UAm50    = round(normrnd(3924, 0.02*3924,[50,1]),0);
Tset50   = round(normrnd(22.5, 0.03*22.5,[50,1]),0);
TDB50    = round(normrnd(1.5,  0.05*1.5, [50,1]),0);

% Helper to simulate one unit for 168 hours under a status policy
simulate_unit = @(policy) ...
    local_sim_unit_minute( tempC_min, policy, ...
                            Prated, Timestep, Tset, TDB, ...
                            UA, UAmass, Vair, DensityAir, Cp_air, Cv_air, C_funiture, V_furniture );

% CASE 1 policy: force off 15–18h
policy_case1 = @(hour_of_day, x1, prev_status, Tset_eff, TDB_eff) deal( ...
                    (hour_of_day>=15 && hour_of_day<18) * 0 + ...
                    (hour_of_day<15 || hour_of_day>=18) * local_hyst(x1(1), prev_status, Tset_eff, TDB_eff) );

% CASE 2 policy: Tset = 18 during 15–18h, otherwise Tset
policy_case2 = @(hour_of_day, x1, prev_status, Tset_eff, TDB_eff) ...
                    local_hyst(x1(1), prev_status, Tset_eff, TDB_eff);

% Build neighborhood loads
for j = 1:50
    % recompute per-house thermal params
    UA = UA50(j); UAmass = UAm50(j);
    R1 = 1/UA; R2 = 1/UAmass;
    Ca = Vair50(j)*Density50(j)*Cp50(j);
    Cm = Vf50(j) * Cf50(j);
    A = [-(1/(R2*Ca)+1/(R1*Ca))  1/(R2*Ca); 1/(R2*Cm)  -1/(R2*Cm)];
    B = [1/(R1*Ca) 1/Ca ;0 0];

    % --- Basecase
    [~, Qb_min] = local_sim_unit_minute_custom(tempC_min, Prated50(j), Timestep, Tset50(j), TDB50(j), A, B, false);
    ngbhd_load(:,1) = ngbhd_load(:,1) + Qb_min;

    % --- Case 1 (force off 15–18h)
    [~, Q1_min] = local_sim_unit_minute_custom(tempC_min, Prated50(j), Timestep, Tset50(j), TDB50(j), A, B, 1);
    ngbhd_load(:,2) = ngbhd_load(:,2) + Q1_min;

    % --- Case 2 (Tset=18 @ 15–18h)
    [~, Q2_min] = local_sim_unit_minute_custom(tempC_min, Prated50(j), Timestep, Tset50(j), TDB50(j), A, B, 2);
    ngbhd_load(:,3) = ngbhd_load(:,3) + Q2_min;
end

% Example plot: day 1 per-minute neighborhood load per case
color = ['k' 'y' 'c' 'b' 'r' 'm' 'g'];
for caseId = 1:3
    Time_hr = (0:1439)/60;
    day_start = 1;
    day_end   = 1440;
    figure(81+caseId)
    set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
    plot(Time_hr, ngbhd_load(day_start:day_end, caseId), color(caseId))
    xlabel('Time (hour)'); ylabel('Load (arb. units)')
    sgtitle(['Neighborhood per-minute load (Day 1) - CASE: ' num2str(caseId)])
end

%% --------- Average Cost Plot for 50 units, Case 1 ----------
load_per_hour = sum(reshape(ngbhd_load(:,2), 60, 168), 1)/60;   % 1x168
cost_per_hour = (load_per_hour .* DA_record(:,10).')/1e6/50;
Time_hr = 0:167;
figure(90)
set(gcf,'DefaultAxesFontSize',10); set(gcf,'DefaultTextFontSize',10)
plot(Time_hr, cost_per_hour, color(1))
xlabel('Time (hour)'); ylabel('Cost ($)')
sgtitle('Average (50 units) hourly load cost over the week (Case 1)')
hold on
for i=1:7
    xline(15 + (i-1)*24, 'r--', 'LineWidth', 2);
    xline(18 + (i-1)*24, 'r--', 'LineWidth', 2);
end

%% --------- Save all figures to /results ----------
outputFolder = fullfile(pwd, 'results');
if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end

figHandles = findall(0, 'Type', 'figure');
for i = 1:length(figHandles)
    figure(figHandles(i));
    filename = sprintf('figure%d.png', i);
    saveas(figHandles(i), fullfile(outputFolder, filename));
end

%% ======== Local helpers ========

function s = local_hyst(Troom, prev_status, Tset, TDB)
% simple on/off with deadband
    if (Troom > (Tset + TDB/2)) && (prev_status == 1)
        s = 0;
    elseif (Troom < (Tset - TDB/2)) && (prev_status == 0)
        s = 1;
    else
        s = prev_status;
    end
end

function [T_min, Q_min] = local_sim_unit_minute_custom(tempC_min, Prated, Timestep, Tset, TDB, A, B, mode)
% mode: 0 base, 1 case1(off 15-18), 2 case2(Tset=18 at 15-18)
    T_min = zeros(10080,1);
    Q_min = zeros(10080,1);
    x0 = [22;22];

    for hh = 1:168
        hour_of_day = mod(hh-1,24);
        switch mode
            case 0
                Tset_eff = Tset;
                force_off = false;
            case 1
                Tset_eff = Tset;
                force_off = (hour_of_day>=15 && hour_of_day<18);
            otherwise
                Tset_eff = (hour_of_day>=15 && hour_of_day<18) * 18 + ...
                           (hour_of_day<15 || hour_of_day>=18) * Tset;
                force_off = false;
        end

        Tout_hour = tempC_min((hh-1)*60 + (1:60));

        Q = zeros(1,60); status = zeros(1,60);
        if     x0(1) > (Tset_eff + TDB/2), status(1)=0;
        elseif x0(1) < (Tset_eff - TDB/2), status(1)=1;
        else,  status(1)=0;
        end
        if force_off, status(1)=0; end

        Q0 = (status(1)==1) * (Prated*Timestep);
        u0 = [Tout_hour(1); Q0];
        Troom_record = zeros(1,60);

        for i=2:60
            DT = A*x0 + B*u0; x1 = x0 + DT;
            if force_off
                status(i) = 0;
            else
                status(i) = local_hyst(x1(1), status(i-1), Tset_eff, TDB);
            end
            if status(i)==1, Q(i)=Prated*Timestep; else, Q(i)=0; end
            x0 = x1; u0 = [Tout_hour(i); Q(i)];
            Troom_record(i-1) = x1(1);
        end

        T_min((hh-1)*60 + (1:60)) = Troom_record(:);
        Q_min((hh-1)*60 + (1:60)) = Q(:);
    end
end
