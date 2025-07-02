%% Tongue direction 100 ms both fd and dr near min protrusion (get similar kinematics)
%% Load data (drinking)
load('20181210_Kinematics_Before.mat');
%% Get tongue tip position
% Min tongue protrusion to spout
ntrials = length(Kinematics.Cranium.points);
Neural = Kinematics.NeuralIndex;

ttptcols = find(contains(Kinematics.ColumnNames.points,'AnteriorM_'));

ttpt = {};
neural = {};
cycles = {}; Cycles = {};
ttAP = {};
Mins = {}; mins = {};
nAP = {}; traj = {}; times = {};
for i = 1:ntrials
    ttpt{i} = Kinematics.Cranium.points{i}(:,ttptcols);
    neural{i} = Neural{i}(:,3);
    npeaks = length(ttpt{i});
    stretch = ttpt{i}(:,1);
    [~,minstretch] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch] = findpeaks(stretch,'MinPeakProminence',10);
    if length(minstretch) > length(maxstretch)
        minstretch = minstretch(1:end-1);
    end
    Mins{i} = ttpt{i}(minstretch,1);
    mins{i} = minstretch;
    for l = 1:length(minstretch)
        ttAP{i}{l} = ttpt{i}(minstretch(l):minstretch(l)+20,1:3);
        nAP{i}{l} = neural{i}(minstretch(l):minstretch(l)+20);
    end
end
traj = horzcat(ttAP{:}); traj = traj(~cellfun('isempty',traj));
times = horzcat(nAP{:}); times = times(~cellfun('isempty',times));
pos1 = find(~cellfun(@isempty,traj));

% Get 3D direction
npz=[0 0 1];
a3 = [];
for t = 1:length(traj)
    v1 = [traj{t}(1,1) traj{t}(1,3) traj{t}(1,2)];
    v2 = [traj{t}(end,1) traj{t}(end,3) traj{t}(end,2)];
    a3(t) = vecangle360(v1,v2,npz);
end
a3=a3'; A3=a3(:);
A3 = A3(pos1);
ControlDir = A3;

dir = ControlDir(ControlDir > -5 & ControlDir < 5); Con = dir;
diridx = find(ControlDir > -5 & ControlDir < 5); ConI = diridx;

% %% Stable Units
% % load RyDataStats2
% inx = scoreW(:,2)>0.99;
% % sum(inx);
% stableUnits = units(inx,:);
% stableM1U = stableUnits(1:34,:); % 1-96

%% Load neural data
M1U = load('20181210b_M1U_sortedspikes.mat');
M1U = struct2cell(M1U);
Units = M1U(Ind_Dr);
%% Get spike counts
Counts_Dr = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(traj)
        if isnan(times{t}(end))
            continue
        elseif isnan(times{t}(1))
            continue
        end
        % 1ms bins for FA %
        edges = times{t}(1) : 0.001 : (times{t}(end) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        Counts_Dr{t}(u,:) = allCounts;
    end
end

%% Load data (feeding)
load('20190228_Kinematics.mat')
%% Get tongue tip position
ttcols = find(contains(Kinematics.ColumnNames.points,'AnteriorM_'));
tonguetip = Kinematics.Cranium.points{2}(:,ttcols);
temptimes = Kinematics.NeuralIndex{2}(:,3);
trial = Kinematics.GapeCycleInfo;
%% Get tongue tip position
% Min tongue protrusion to spout
ntrials = length(Kinematics.Cranium.points);
Neural = Kinematics.NeuralIndex;

ttptcols = find(contains(Kinematics.ColumnNames.points,'AnteriorM_'));

ttpt = {};
neural = {};
cycles = {}; Cycles = {};
ttAP = {};
Mins = {}; mins = {};
nAP = {}; traj = {}; times = {};
for i = 1:ntrials
    ttpt{i} = Kinematics.Cranium.points{i}(:,ttptcols);
    neural{i} = Neural{i}(:,3);
    npeaks = length(ttpt{i});
    stretch = ttpt{i}(:,1);
    [~,minstretch] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch] = findpeaks(stretch,'MinPeakProminence',10);
    if length(minstretch) > length(maxstretch)
        minstretch = minstretch(1:end-1);
    end
    Mins{i} = ttpt{i}(minstretch,1);
    mins{i} = minstretch;
    for l = 1:length(minstretch)
        ttAP{i}{l} = ttpt{i}(minstretch(l):minstretch(l)+20,1:3);
        nAP{i}{l} = neural{i}(minstretch(l):minstretch(l)+20);
    end
end
traj = horzcat(ttAP{:}); traj = traj(~cellfun('isempty',traj));
times = horzcat(nAP{:}); times = times(~cellfun('isempty',times));
pos1 = find(~cellfun(@isempty,traj));

% Get 3D direction
npz=[0 0 1];
a3 = [];
for t = 1:length(traj)
    v1 = [traj{t}(1,1) traj{t}(1,3) traj{t}(1,2)];
    v2 = [traj{t}(end,1) traj{t}(end,3) traj{t}(end,2)];
    a3(t) = vecangle360(v1,v2,npz);
end
a3=a3'; A3=a3(:);
A3 = A3(pos1);
ControlDir = A3;

dir = ControlDir(ControlDir > -5 & ControlDir < 5); Con = dir;
diridx = find(ControlDir > -5 & ControlDir < 5); ConI = diridx;
%% Load neural data
M1U = load('20190228_M1U_sortedspikes.mat');
M1U = struct2cell(M1U);
Units = M1U(Ind_Fd);
%% Get spike counts
Counts_Fd = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(traj)
        if isnan(times{t}(end))
            continue
        elseif isnan(times{t}(1))
            continue
        end
        % 1ms bins for FA %
        edges = times{t}(1) : 0.001 : (times{t}(end) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        Counts_Fd{t}(u,:) = allCounts;
    end
end

%% Combine drinking and feeding
Counts_Fd = Counts_Fd(~cellfun('isempty',Counts_Fd)) ; 
Counts_Dr = Counts_Dr(~cellfun('isempty',Counts_Dr)) ; 

F = cell(1,length(Counts_Fd)); F(:) = {'Feed'};
D = cell(1,length(Counts_Dr)); D(:) = {'Drink'};
vecI = [F D];
condition = vecI'; %% directions

data = [Counts_Fd Counts_Dr];
Data = struct();
for i = 1:length(data)
    Data(i).data = data{i};
    Data(i).condition = condition{i};
%     {data{i}; vecI{i}};
    Data(i).type = 'traj';
    Data(i).epochStarts = 1;
    Data(i).epochColors = [0.598590500684565,0.668487146264776,0.894564090918275];
end
%% Dimensionality reduction
DataHigh(Data, 'DimReduce');
