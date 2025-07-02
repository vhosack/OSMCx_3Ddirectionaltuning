%% Jaw kinematics during Drinking
%% Load Data
Before = load('20181210_Kinematics_Before.mat');%%%
After = load('20181210_Kinematics_After.mat');
% copy and paste paired dataset paths
%% Min tongue protrusion to spout
ntrials = length(Before.Kinematics.Cranium.points);
Neural = Before.Kinematics.NeuralIndex;
JawTranslation = Before.Kinematics.TMJ;

ttptcols = find(contains(Before.Kinematics.ColumnNames.points,'AnteriorM_'));

ttpt = {};
neural = {};
cycles = {}; Cycles = {}; times = {};
ttAP = {}; JawKin = {};
Mins = {}; mins = {};
for i = 1:ntrials
    ttpt{i} = Before.Kinematics.Cranium.points{i}(:,ttptcols);
    neural{i} = Neural{i}(:,3);
    npeaks = length(ttpt{i});
    stretch = ttpt{i}(:,1);
    [~,minstretch] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch] = findpeaks(stretch,'MinPeakProminence',10);
    if length(minstretch) > length(maxstretch)
        minstretch = minstretch(1:end-1);
    end
%     Mins = [ttpt{i}(:,1), ttpt{i}(:,2), ttpt{i}(:,3)]; mins = Mins(minstretch,:);
%     Maxs = [ttpt{i}(:,1), ttpt{i}(:,2), ttpt{i}(:,3)]; maxs = Maxs(maxstretch,:);
    Mins{i} = ttpt{i}(minstretch,1);
    mins{i} = minstretch;
    ttAP{i} = ttpt{i}(minstretch:maxstretch,1);
    JawKin{i} = JawTranslation{i}(minstretch:maxstretch,6);
    if isempty(JawKin{i})
        continue
    end
    ttAP{i} = resample(ttAP{i},100,length(ttAP{i}));
    JawKin{i} = resample(JawKin{i},100,length(JawKin{i}));
end
ttAP = ttAP(~cellfun('isempty',ttAP));
allttAP = cat(3,ttAP{:});
meanttAP = mean(allttAP, 3);
JawKin = JawKin(~cellfun('isempty',JawKin));
allJawKin = cat(3,JawKin{:});
meanJawKin = mean(allJawKin, 3);
stdJawKin = std(allJawKin, [], 3);
lower = meanJawKin - stdJawKin;
upper = meanJawKin + stdJawKin;
%% PLot Jaw Pitch and TT Ant-Post (1 trial)
t = 0.005:0.005:11.005;
figure;
Pitch = JawTranslation{1}(:,6);
plot(t,Pitch,'k'); ylabel('Jaw Pitch'); title('Drinking'); ylim([-25 0]); xlabel('Time (s)'); xlim([1 10]); hold on;
yyaxis right
TT = ttpt{1}(:,1);
plot(t,TT); ylabel('Anterior-Posterior position'); ylim([20 65]);
scatter(mins{1}/200,Mins{1},'^','filled');

%% Plot Jaw Pitch and TT Ant-Post (mean)
figure;
t = 1:100;
plot(meanJawKin); ylabel('Jaw Pitch'); hold on;
yyaxis right
plot(meanttAP); ylabel('Ant-Post TT pos');

%% Plot all variables (mean)
figure;
t = 1:100;
subplot(3,1,1); plot(meanJawKin(:,1)); ylabel('x-translation'); hold on;
% edges = [t, fliplr(t)]; sem = [upper(:,1)', fliplr(lower(:,1))']; fill(edges,sem,'b');
subplot(3,1,2); plot(meanJawKin(:,2)); ylabel('y-translation'); hold on;
% edges = [t, fliplr(t)]; sem = [upper(:,2)', fliplr(lower(:,2))']; fill(edges,sem,'b');
subplot(3,1,3); plot(meanJawKin(:,3)); ylabel('z-translation'); xlabel('% of Cycle');  hold on;
% edges = [t, fliplr(t)]; sem = [upper(:,3)', fliplr(lower(:,3))']; fill(edges,sem,'b');

%% Feeding
%% Load Data
Control = load('20190228_Kinematics.mat');
NerveBlock = load('20190227_Kinematics.mat');
Kinematics = [Control NerveBlock];
% copy and paste paired dataset paths
%% Jaw Pitch vs Ant-Post TT pos
t = 0.005:0.005:11.005;
figure;
JawFeed = Control.Kinematics.TMJ;
Pitch = JawFeed{1}(:,6);
plot(t,Pitch,'k'); ylabel('Jaw Pitch'); title('Feeding'); xlabel('Time (s)'); xlim([1 10]); hold on;
yyaxis right
ttptcols_Fd = find(contains(Control.Kinematics.ColumnNames.points,'AnteriorM_'));
ttpt_Fd = Control.Kinematics.Cranium.points{1}(:,ttptcols_Fd);
TT = ttpt_Fd(:,1);
plot(t,TT); ylabel('Anterior-Posterior position'); %ylim([23 65]);

% trial = Control.Kinematics.TrialNames{1};
% tempdata = Control.Kinematics.GapeCycleInfo(strcmp(trial,Control.Kinematics.GapeCycleInfo.Trialname),:);
% mims = tempdata.MinGape;
% Mims = TT(mims);
% scatter(mims/200,Mims,'^','filled');
