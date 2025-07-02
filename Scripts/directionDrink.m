%% Get directions and spike data (drinking)
%% Load Data
Before = load('20181210_Kinematics_Before.mat');
After = load('20181210_Kinematics_After.mat');
% copy and paste paired dataset paths
%% Comparing Directions (Before)
right = find(Before.Kinematics.SpoutNumber == 2);
middle = find(Before.Kinematics.SpoutNumber == 1);
left = find(Before.Kinematics.SpoutNumber == 3);

Right = Before.Kinematics.Cranium.points(right);
Middle = Before.Kinematics.Cranium.points(middle);
Left = Before.Kinematics.Cranium.points(left);

ttptcols = find(contains(Before.Kinematics.ColumnNames.points,'AnteriorM_'));

figure;
% Right trials
ntrials = length(Right);
Neural_r = Before.Kinematics.NeuralIndex(right);
ttpt_r = {};
neural_r = {};
cycles_r = {}; Cycles_r = {};
times_r = {};
allintrial_r = {}; numtrials_r=[]; avglength_r = {};
for i = 1:ntrials
    ttpt_r{i} = Right{i}(:,ttptcols);
    neural_r{i} = Neural_r{i}(:,3);
    npeaks = length(ttpt_r{i});
    stretch = ttpt_r{i}(:,2);
    [~,minstretch_r] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_r] = findpeaks(stretch,'MinPeakProminence',10);
    if length(minstretch_r) == length(maxstretch_r)
        allintrial_r{i} = abs(minstretch_r - maxstretch_r);
        numtrials_r(i) = length(allintrial_r{i});
        avglength_r{i} = allintrial_r{i}(allintrial_r{i}<200);
    elseif length(minstretch_r) > length(maxstretch_r)
        minstretch = minstretch_r(1:end-1);
        allintrial_r{i} = abs(minstretch - maxstretch_r);
        numtrials_r(i) = length(allintrial_r{i});
        avglength_r{i} = allintrial_r{i}(allintrial_r{i}<200);
    end
    Mins = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    mins = Mins(minstretch_r,:);
    Maxs = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    maxs = Maxs(maxstretch_r,:);
    minlength = length(maxstretch_r);
    for k = 1:minlength
        origin = maxstretch_r(k);
        cycles_r{i,k}(:,1) = ttpt_r{i}(origin-100:origin+100,1);
        cycles_r{i,k}(:,2) = ttpt_r{i}(origin-100:origin+100,2);
        cycles_r{i,k}(:,3) = ttpt_r{i}(origin-100:origin+100,3);
%         Visualize trajectories
        plot3(cycles_r{i,k}(:,1), cycles_r{i,k}(:,2), cycles_r{i,k}(:,3), 'y');
        hold on;
        times_r{i,k} = neural_r{i}(maxstretch_r(k));
    end
    Cycles_r = cycles_r(~cellfun(@isempty,cycles_r))';
    Times_r = times_r(~cellfun(@isempty,times_r))';
end

% Middle trials
ntrials = length(Middle);
Neural_m = Before.Kinematics.NeuralIndex(middle);
ttpt_m = {};
neural_m = {};
cycles_m = {}; Cycles_m = {};
times_m = {};
numtrials_m=[]; avglength_m = {};
for i = 1:ntrials
    ttpt_m{i} = Middle{i}(:,ttptcols);
    neural_m{i} = Neural_m{i}(:,3);
    npeaks = length(ttpt_m{i});
    stretch = ttpt_m{i}(:,2);
    [~,minstretch_m] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_m] = findpeaks(stretch,'MinPeakProminence',10);
    if length(minstretch_m) == length(maxstretch_m)
        allintrial_m{i} = abs(minstretch_m - maxstretch_m);
        numtrials_m(i) = length(allintrial_m{i});
        avglength_m{i} = allintrial_m{i}(allintrial_m{i}<200);
    elseif length(minstretch_m) > length(maxstretch_m)
        minstretch = minstretch_m(1:end-1);
        allintrial_m{i} = abs(minstretch - maxstretch_m);
        numtrials_m(i) = length(allintrial_m{i});
        avglength_m{i} = allintrial_m{i}(allintrial_m{i}<200);
    end
    npeaks = length(minstretch_m);
    Mins = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    maxs = Maxs(maxstretch_m,:);
    minlength = length(maxstretch_m);
    for k = 1:minlength
        origin = maxstretch_m(k);
        cycles_m{i,k}(:,1) = ttpt_m{i}(origin-100:origin+100,1);
        cycles_m{i,k}(:,2) = ttpt_m{i}(origin-100:origin+100,2);
        cycles_m{i,k}(:,3) = ttpt_m{i}(origin-100:origin+100,3);
%         Visualize trajectories
        plot3(cycles_m{i,k}(:,1), cycles_m{i,k}(:,2), cycles_m{i,k}(:,3), 'b');
        hold on;
        times_m{i,k} = neural_m{i}(maxstretch_m(k));
    end
    Cycles_m = cycles_m(~cellfun(@isempty,cycles_m))';
    Times_m = times_m(~cellfun(@isempty,times_m))';
end

% Left trials
ntrials = length(Left);
Neural_l = Before.Kinematics.NeuralIndex(left);
ttpt_l = {};
neural_l = {};
cycles_l = {}; Cycles_l = {};
times_l = {};
numtrials_l = []; avglength_l = {};
for i = 1:ntrials
    ttpt_l{i} = Left{i}(:,ttptcols);
    neural_l{i} = Neural_l{i}(:,3);
    npeaks = length(ttpt_l{i});
    stretch = ttpt_l{i}(:,2);
    [~,minstretch_l] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_l] = findpeaks(stretch,'MinPeakProminence',10);
    if length(minstretch_l) == length(maxstretch_l)
        allintrial_l{i} = abs(minstretch_l - maxstretch_l);
        numtrials_l(i) = length(allintrial_l{i});
        avglength_l{i} = allintrial_l{i}(allintrial_l{i}<200);
    elseif length(minstretch_l) > length(maxstretch_l)
        minstretch = minstretch_l(1:end-1);
        allintrial_l{i} = abs(minstretch - maxstretch_l);
        numtrials_l(i) = length(allintrial_l{i});
        avglength_l{i} = allintrial_l{i}(allintrial_l{i}<200);
    end
    npeaks = length(minstretch_l);
    Mins = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    maxs = Maxs(maxstretch_l,:);
    minlength = length(maxstretch_l);
    for k = 1:minlength
        origin = maxstretch_l(k);
        cycles_l{i,k}(:,1) = ttpt_l{i}(origin-100:origin+100,1);
        cycles_l{i,k}(:,2) = ttpt_l{i}(origin-100:origin+100,2);
        cycles_l{i,k}(:,3) = ttpt_l{i}(origin-100:origin+100,3);
%         Visualize trajectories
        plot3(cycles_l{i,k}(:,1), cycles_l{i,k}(:,2), cycles_l{i,k}(:,3), 'g');
        hold on;
        times_l{i,k} = neural_l{i}(maxstretch_l(k));
    end
    Cycles_l = cycles_l(~cellfun(@isempty,cycles_l))';
    Times_l = times_l(~cellfun(@isempty,times_l))';
end
IndBefore = [Times_l'; Times_m'; Times_r'];
CyclesBefore = {Cycles_l;Cycles_m;Cycles_r};
NeuralBefore = {Times_l;Times_m;Times_r};
%% Neural (Before) - run each area separately
%% MIo
M1F = load('20190621b_M1F_sortedspikes.mat');
M1U = load('20190621b_M1U_sortedspikes.mat');
M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
Units = [M1F; M1U];
% Units = M1U; %% run individual area

spikecount = {};
meanfiring = {};
timespikes = {};
Counts = {};
timeframe = {};
spiketimes = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for dir = 1:length(NeuralBefore)
        trials = NeuralBefore{dir};
        trials = trials(~cellfun(@isempty,trials));
        cycles = cell2mat(trials);
        Cycles{u}{dir} = cycles;
        first = cycles-0.25;  %% time in s
        last = cycles+0.25;   %%
        timeframe{dir} = [first' last'];
        for c = 1:length(cycles)
            
            edges = [];
            % 1ms bins for FA %
            edges = first(c) : 0.001 : (last(c) + 0.001); %
            allCounts = zeros(1, length(edges) - 1); %
            [counts, values] = histcounts(temp, edges); %
            allCounts = allCounts + counts; %
            Counts{dir}{c}(u,:) = allCounts;

            inx = find(temp >= first(c) & temp <= last(c));
            if inx == 0
                spikecount{dir}(c,u) = [];
            end
            spikecount{dir}(c,u) = length(inx);
            Ts = []; Ts = temp(inx);
            if ~isempty(Ts)
                Ts = Ts - temp(1);
            end
            timespikes{u}{dir}{c} = Ts;
            spiketimes{dir}{c,u} = Ts;
            time = last(c) - first(c);
            meanfiring{dir}(c,u) = spikecount{dir}(c,u)/time;
        end
    end
end

for d = 1:3
    cycles = meanfiring{d};
    ncycles = length(cycles);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
    end
end
M1B = meanfr;

% Tuning Curves
% yFit = [];
% yFit3 = [];
% angle = 1:3; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = M1B(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
% %     yFit(:,unit) = cos_fun(p,-30:30);
%     yFit3(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit3(:,unit));
% end

% for u = 1:length(Units)
%     figure;
%     plot(angle,M1FB(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
%     %plot([-30:10:30],yFit8(:,u),'r-.') % cosine fit
%     %legend('Actual Data','Cosine Fit')
% end
%% SIo
S1F = load('20190621b_S1F_sortedspikes.mat');
S1U = load('20190621b_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = [S1F; S1U];
% Units = S1U; %% run individual area

spikecount = {};
meanfiring = {};
timespikes = {};
Counts = {};
timeframe = {};
spiketimes = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for dir = 1:length(NeuralBefore)
        trials = NeuralBefore{dir};
        trials = trials(~cellfun(@isempty,trials));
        cycles = cell2mat(trials);
        Cycles{u}{dir} = cycles;
        first = cycles-0.25;  %% time in s
        last = cycles+0.25;   %%
        timeframe{dir} = [first' last'];
        for c = 1:length(cycles)

            % 1ms bins for FA %
            edges = first(c) : 0.001 : (last(c) + 0.001); %
            allCounts = zeros(1, length(edges) - 1); %
            [counts, values] = histcounts(temp, edges); %
            allCounts = allCounts + counts; %
            Counts{dir}{c}(u,:) = allCounts;

            inx = find(temp >= first(c) & temp <= last(c));
            if inx == 0
                spikecount{dir}(c,u) = [];
            end
            spikecount{dir}(c,u) = length(inx);
            Ts = []; Ts = temp(inx);
            if ~isempty(Ts)
                Ts = Ts - temp(1);
            end
            timespikes{u}{dir}{c} = Ts;
            spiketimes{dir}{c,u} = Ts;
            time = last(c) - first(c);
            meanfiring{dir}(c,u) = spikecount{dir}(c,u)/time;
        end
    end
end

for d = 1:3
    cycles = meanfiring{d};
    ncycles = length(cycles);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
    end
end
S1B = meanfr;

% % Tuning Curves
% yFit = [];
% yFit3 = [];
% angle = 1:3; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = S1B(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
% %     yFit(:,unit) = cos_fun(p,-30:30);
%     yFit3(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit3(:,unit));
% end

% for u = 1:length(Units)
%     figure;
%     plot(angle,S1B(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
%     %plot([-30:10:30],yFit8(:,u),'r-.') % cosine fit
%     %legend('Actual Data','Cosine Fit')
% end
%% Comparing Directions (After)
right = find(After.Kinematics.SpoutNumber == 2);
middle = find(After.Kinematics.SpoutNumber == 1);
left = find(After.Kinematics.SpoutNumber == 3);

Right = After.Kinematics.Cranium.points(right);
Middle = After.Kinematics.Cranium.points(middle);
Left = After.Kinematics.Cranium.points(left);

ttptcols = find(contains(After.Kinematics.ColumnNames.points,'AnteriorM_'));

% Right trials
ntrials = length(Right);
Neural_r = After.Kinematics.NeuralIndex(right);
ttpt_r = {};
neural_r = {};
cycles_r = {};
Cycles_r = {};
times_r = {};
for i = 1:ntrials
    ttpt_r{i} = Right{i}(:,ttptcols);
    neural_r{i} = Neural_r{i}(:,3);
    npeaks = length(ttpt_r{i});
    stretch = ttpt_r{i}(:,2);
    [~,minstretch_r] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_r] = findpeaks(stretch,'MinPeakProminence',10);
    Mins = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    mins = Mins(minstretch_r,:);
    Maxs = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    maxs = Maxs(maxstretch_r,:);
    minlength = length(maxstretch_r);
    for k = 1:minlength
        origin = maxstretch_r(k);
        cycles_r{i,k}(:,1) = ttpt_r{i}(origin-100:origin+100,1);
        cycles_r{i,k}(:,2) = ttpt_r{i}(origin-100:origin+100,2);
        cycles_r{i,k}(:,3) = ttpt_r{i}(origin-100:origin+100,3);
%         % Visualize trajectories
%         plot3(cycles_r{i,k}(:,1), cycles_r{i,k}(:,2), cycles_r{i,k}(:,3), 'y');
%         hold on;
        times_r{i,k} = neural_r{i}(maxstretch_r(k));
    end
    Cycles_r = cycles_r(~cellfun(@isempty,cycles_r))';
    Times_r = times_r(~cellfun(@isempty,times_r))';
end

% Middle trials
ntrials = length(Middle);
Neural_m = After.Kinematics.NeuralIndex(middle);
ttpt_m = {};
neural_m = {};
cycles_m = {};
Cycles_m = {};
times_m = {};
for i = 1:ntrials
    ttpt_m{i} = Middle{i}(:,ttptcols);
    neural_m{i} = Neural_m{i}(:,3);
    npeaks = length(ttpt_m{i});
    stretch = ttpt_m{i}(:,2);
    [~,minstretch_m] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_m] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_m);
    Mins = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    maxs = Maxs(maxstretch_m,:);
    minlength = length(maxstretch_m);
    for k = 1:minlength
        origin = maxstretch_m(k);
        cycles_m{i,k}(:,1) = ttpt_m{i}(origin-100:origin+100,1);
        cycles_m{i,k}(:,2) = ttpt_m{i}(origin-100:origin+100,2);
        cycles_m{i,k}(:,3) = ttpt_m{i}(origin-100:origin+100,3);
%         % Visualize trajectories
%         plot3(cycles_m{i,k}(:,1), cycles_m{i,k}(:,2), cycles_m{i,k}(:,3), 'b');
%         hold on;
        times_m{i,k} = neural_m{i}(maxstretch_m(k));
    end
    Cycles_m = cycles_m(~cellfun(@isempty,cycles_m))';
    Times_m = times_m(~cellfun(@isempty,times_m))';
end

% Left trials
ntrials = length(Left);
Neural_l = After.Kinematics.NeuralIndex(left);
ttpt_l = {};
neural_l = {};
cycles_l = {};
Cycles_l = {};
times_l = {};
for i = 1:ntrials
    ttpt_l{i} = Left{i}(:,ttptcols);
    neural_l{i} = Neural_l{i}(:,3);
    npeaks = length(ttpt_l{i});
    stretch = ttpt_l{i}(:,2);
    [~,minstretch_l] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_l] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_l);
    Mins = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    maxs = Maxs(maxstretch_l,:);
    minlength = length(maxstretch_l);
    for k = 1:minlength
        origin = maxstretch_l(k);
        cycles_l{i,k}(:,1) = ttpt_l{i}(origin-100:origin+100,1);
        cycles_l{i,k}(:,2) = ttpt_l{i}(origin-100:origin+100,2);
        cycles_l{i,k}(:,3) = ttpt_l{i}(origin-100:origin+100,3);
%         % Visualize trajectories
%         plot3(cycles_r{i,k}(:,1), cycles_r{i,k}(:,2), cycles_r{i,k}(:,3), 'g');
%         hold on;
        times_l{i,k} = neural_l{i}(maxstretch_l(k));
    end
    Cycles_l = cycles_l(~cellfun(@isempty,cycles_l))';
    Times_l = times_l(~cellfun(@isempty,times_l))';
end
IndAfter = [Times_l'; Times_m'; Times_r'];
CyclesAfter = {Cycles_l;Cycles_m;Cycles_r};
NeuralAfter = {Times_l;Times_m;Times_r};
%% Neural (After) - run each area separately
%% MIo
M1F = load('20190621a_M1F_sortedspikes.mat');
M1U = load('20190621a_M1U_sortedspikes.mat');
M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
Units = [M1F; M1U];
% Units = M1U; %% run individual area

spikecount = {};
meanfiring = {};
timespikes = {};
Counts = {};
timeframe = {};
spiketimes = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for dir = 1:length(NeuralAfter)
        trials = NeuralAfter{dir};
        trials = trials(~cellfun(@isempty,trials));
        cycles = cell2mat(trials);
        Cycles{u}{dir} = cycles;
        first = cycles-0.25;  %% time in s
        last = cycles+0.25;   %%
        timeframe{dir} = [first' last'];
        for c = 1:length(cycles)

            % 1ms bins for FA %
            edges = first(c) : 0.001 : (last(c) + 0.001); %
            allCounts = zeros(1, length(edges) - 1); %
            [counts, values] = histcounts(temp, edges); %
            allCounts = allCounts + counts; %
            Counts{dir}{c}(u,:) = allCounts;

            inx = find(temp >= first(c) & temp <= last(c));
            if inx == 0
                spikecount{dir}(c,u) = [];
            end
            spikecount{dir}(c,u) = length(inx);
            Ts = []; Ts = temp(inx);
            if ~isempty(Ts)
                Ts = Ts - temp(1);
            end
            timespikes{u}{dir}{c} = Ts;
            spiketimes{dir}{c,u} = Ts;
            time = last(c) - first(c);
            meanfiring{dir}(c,u) = spikecount{dir}(c,u)/time;
        end
    end
end

for d = 1:3
    cycles = meanfiring{d};
    ncycles = length(cycles);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
    end
end
M1A = meanfr;

% % Tuning Curves
% yFit = [];
% yFit3 = [];
% angle = 1:3; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = M1A(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
% %     yFit(:,unit) = cos_fun(p,-30:30);
%     yFit3(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit3(:,unit));
% end

% for u = 1:length(Units)
%     figure;
%     plot(angle,M1B(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
%     %plot([-30:10:30],yFit8(:,u),'r-.') % cosine fit
%     %legend('Actual Data','Cosine Fit')
% end
%% SIo
S1F = load('20190621a_S1F_sortedspikes.mat');
S1U = load('20190621a_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = [S1F; S1U];
% Units = S1U; %% run individual area

spikecount = {};
meanfiring = {};
timespikes = {};
Counts = {};
timeframe = {};
spiketimes = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for dir = 1:length(NeuralAfter)
        trials = NeuralAfter{dir};
        trials = trials(~cellfun(@isempty,trials));
        cycles = cell2mat(trials);
        Cycles{u}{dir} = cycles;
        first = cycles-0.25;  %% time in s
        last = cycles+0.25;   %%
        timeframe{dir} = [first' last'];
        for c = 1:length(cycles)

            % 1ms bins for FA %
            edges = first(c) : 0.001 : (last(c) + 0.001); %
            allCounts = zeros(1, length(edges) - 1); %
            [counts, values] = histcounts(temp, edges); %
            allCounts = allCounts + counts; %
            Counts{dir}{c}(u,:) = allCounts;

            inx = find(temp >= first(c) & temp <= last(c));
            if inx == 0
                spikecount{dir}(c,u) = [];
            end
            spikecount{dir}(c,u) = length(inx);
            Ts = []; Ts = temp(inx);
            if ~isempty(Ts)
                Ts = Ts - temp(1);
            end
            timespikes{u}{dir}{c} = Ts;
            spiketimes{dir}{c,u} = Ts;
            time = last(c) - first(c);
            meanfiring{dir}(c,u) = spikecount{dir}(c,u)/time;
        end
    end
end

for d = 1:3
    cycles = meanfiring{d};
    ncycles = length(cycles);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
    end
end
S1A = meanfr;

% % Tuning Curves
% yFit = [];
% yFit3 = [];
% angle = 1:3; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = S1A(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
% %     yFit(:,unit) = cos_fun(p,-30:30);
%     yFit3(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit3(:,unit));
% end

% for u = 1:length(Units)
%     figure;
%     plot(angle,S1B(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
%     %plot([-30:10:30],yFit8(:,u),'r-.') % cosine fit
%     %legend('Actual Data','Cosine Fit')
% end
