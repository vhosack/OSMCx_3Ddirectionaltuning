%% Obtain spike times for psth creation
%% Load Data (Drinking)
Before = load('20190621_Kinematics_Before.mat');
After = load('20190621_Kinematics_After.mat');
%% Comparing Directions (Before)
right = find(Before.Kinematics.SpoutNumber == 2);
middle = find(Before.Kinematics.SpoutNumber == 1);
left = find(Before.Kinematics.SpoutNumber == 3);

Right = Before.Kinematics.Cranium.points(right);
Middle = Before.Kinematics.Cranium.points(middle);
Left = Before.Kinematics.Cranium.points(left);

ttptcols = find(contains(Before.Kinematics.ColumnNames.points,'AnteriorM_'));

% Right trials
ntrials = length(Right);
Neural_r = Before.Kinematics.NeuralIndex(right);
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
    %plot3(mins(:,1), mins(:,2), mins(:,3), 'ko');
    Maxs = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    maxs = Maxs(maxstretch_r,:);
    %plot3(maxs(:,1), maxs(:,2), maxs(:,3), 'mo');
    minlength = length(maxstretch_r);
    for k = 1:minlength
        origin = maxstretch_r(k);
        cycles_r{i,k}(:,1) = ttpt_r{i}(origin-100:origin+100,1);
        cycles_r{i,k}(:,2) = ttpt_r{i}(origin-100:origin+100,2);
        cycles_r{i,k}(:,3) = ttpt_r{i}(origin-100:origin+100,3);
        plot3(cycles_r{i,k}(:,1), cycles_r{i,k}(:,2), cycles_r{i,k}(:,3), 'b');
        hold on;
        times_r{i,k} = neural_r{i}(maxstretch_r(k));
    end
end
Cycles_r = cycles_r(~cellfun(@isempty,cycles_r));
Times_r = times_r(~cellfun(@isempty,times_r));
% title('Tongue Tip Trajectory')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% xlim([20 70])
% ylim([-50 10])
% zlim([-30 30])

% Middle trials
ntrials = length(Middle);
Neural_m = Before.Kinematics.NeuralIndex(middle);
ttpt_m = {};
neural_m = {};
cycles_m = {};
Cycles_m = {};
times_m = {};
for i = 1:ntrials
    ttpt_m{i} = Middle{i}(:,ttptcols);
    neural_m{i} = Neural_m{i}(:,3);
    npeaks = length(ttpt_m{i});
    %plot3(ttpt_m{i}(:,3), ttpt_m{i}(:,2), ttpt_m{i}(:,1), 'b');
    stretch = ttpt_m{i}(:,2);
    [~,minstretch_m] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_m] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_m);
    Mins = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    mins = Mins(minstretch_m,:);
    %plot3(mins(:,1), mins(:,2), mins(:,3), 'ko');
    Maxs = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    maxs = Maxs(maxstretch_m,:);
    %plot3(maxs(:,1), maxs(:,2), maxs(:,3), 'mo');
    minlength = length(maxstretch_m);
    for k = 1:minlength
        origin = maxstretch_m(k);
        cycles_m{i,k}(:,1) = ttpt_m{i}(origin-100:origin+100,1);
        cycles_m{i,k}(:,2) = ttpt_m{i}(origin-100:origin+100,2);
        cycles_m{i,k}(:,3) = ttpt_m{i}(origin-100:origin+100,3);
        plot3(cycles_m{i,k}(:,1), cycles_m{i,k}(:,2), cycles_m{i,k}(:,3), 'b');
        times_m{i,k} = neural_m{i}(maxstretch_m(k));
    end
end
Cycles_m = cycles_m(~cellfun(@isempty,cycles_m));
Times_m = times_m(~cellfun(@isempty,times_m));

% Left trials
ntrials = length(Left);
Neural_l = Before.Kinematics.NeuralIndex(left);
ttpt_l = {};
neural_l = {};
cycles_l = {};
Cycles_l = {};
times_l = {};
for i = 1:ntrials
    ttpt_l{i} = Left{i}(:,ttptcols);
    neural_l{i} = Neural_l{i}(:,3);
    npeaks = length(ttpt_l{i});
    %plot3(ttpt_l{i}(:,3), ttpt_l{i}(:,2), ttpt_l{i}(:,1), 'g');
    stretch = ttpt_l{i}(:,2);
    [~,minstretch_l] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_l] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_l);
    Mins = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    mins = Mins(minstretch_m,:);
    %plot3(mins(:,1), mins(:,2), mins(:,3), 'ko');
    Maxs = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    maxs = Maxs(maxstretch_l,:);
    %plot3(maxs(:,1), maxs(:,2), maxs(:,3), 'mo');
    minlength = length(maxstretch_l);
    for k = 1:minlength
        origin = maxstretch_l(k);
        cycles_l{i,k}(:,1) = ttpt_l{i}(origin-100:origin+100,1);
        cycles_l{i,k}(:,2) = ttpt_l{i}(origin-100:origin+100,2);
        cycles_l{i,k}(:,3) = ttpt_l{i}(origin-100:origin+100,3);
        plot3(cycles_l{i,k}(:,1), cycles_l{i,k}(:,2), cycles_l{i,k}(:,3), 'b');
        times_l{i,k} = neural_l{i}(maxstretch_l(k));
    end
end
Cycles_l = cycles_l(~cellfun(@isempty,cycles_l));
Times_l = times_l(~cellfun(@isempty,times_l));

% figure;
% plot3(begpoints_r(:,1),begpoints_r(:,2),begpoints_r(:,3),'ro');
% title('Endpoints Before')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% hold on;
% plot3(begpoints_m(:,1),begpoints_m(:,2),begpoints_m(:,3),'bo');
% plot3(begpoints_l(:,1),begpoints_l(:,2),begpoints_l(:,3),'go');

CyclesBefore = {Cycles_l;Cycles_m;Cycles_r};
NeuralBefore = {Times_l;Times_m;Times_r};
%% Neural (Before)
%% Motor Cortex
% % Select dataset
% M1F = load('20181210b_M1F_sortedspikes.mat');
% M1U = load('20181210b_M1U_sortedspikes.mat');
% M1F = struct2cell(M1F); 
% M1U = struct2cell(M1U);
% Units = {[M1F; M1U]};
% Units = Units{1};
% % Units = M1U;

% S1F = load('20181210b_S1F_sortedspikes.mat');
% S1U = load('20181210b_S1U_sortedspikes.mat');
% S1F = struct2cell(S1F); 
% S1U = struct2cell(S1U);
% Units = {[S1F; S1U]};
% Units = Units{1};
% % Units = S1U;

% M1F = load('20190621b_M1F_sortedspikes.mat');
% M1U = load('20190621b_M1U_sortedspikes.mat');
% M1F = struct2cell(M1F); 
% M1U = struct2cell(M1U);
% Units = {[M1F; M1U]};
% Units = Units{1};
% % Units = M1U;

S1F = load('20190621b_S1F_sortedspikes.mat');
S1U = load('20190621b_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = {[S1F; S1U]};
Units = Units{1};
% Units = S1U;

timespikes = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for dir = 1:length(NeuralBefore)
        cycles = NeuralBefore{dir};
        cycles = cell2mat(cycles);
        Cycles{u}{dir} = cycles;
        first = cycles-0.5;
        last = cycles+0.5;
        for c = 1:length(cycles)
            inx = find(temp >= first(c) & temp <= last(c));
            Ts = []; Ts = temp(inx);
            Ts = Ts - cycles(c);
            timespikes{u}{dir}{c} = Ts;
        end
    end
end
%% Motor Cortex
% M1F = load('20181210b_M1F_sortedspikes.mat');
% M1U = load('20181210b_M1U_sortedspikes.mat');
% M1F = struct2cell(M1F); 
% M1U = struct2cell(M1U);
% Units = {[M1F; M1U]};
% Units = Units{1};
% % Units = M1U;

S1F = load('20181210b_S1F_sortedspikes.mat');
S1U = load('20181210b_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = {[S1F; S1U]};
Units = Units{1};
% Units = S1U;

spikecount = {};
meanfiring = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for dir = 1:length(NeuralBefore)
        cycles = NeuralBefore{dir};
        cycles = cell2mat(cycles);
        Cycles{u}{dir} = cycles;
        first = cycles-0.05;
        last = cycles+0.05;
        for c = 1:length(cycles)
            inx = find(temp >= first(c) & temp <= last(c));
            if inx == 0
                spikecount{dir}(c,u) = [];
            end
            spikecount{dir}(c,u) = length(inx);
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
M1FB = meanfr;

% Tuning Curves
% yFit = [];
% yFit3 = [];
% angle = 1:3; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = M1FB(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
% %     yFit(:,unit) = cos_fun(p,-40:40);
%     yFit3(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit3(:,unit));
% end

% for u = 1:length(Units)
%     figure;
%     plot(angle,M1FB(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
%     %plot([-40:10:40],yFit8(:,u),'r-.') % cosine fit
%     %legend('Actual Data','Cosine Fit')
% end

%Preferred Directions
pref = zeros(length(Units),1);
for unit = 1:length(Units)
    [p, d] = max(M1FB(:,unit));
    pref(unit) = d;
end

figure;
histogram(pref,'normalization','pdf');
xlabel('Preferred Direction');
ylabel('Proportion of Neurons');
title('Histogram of Preferred Directions of Motor Neurons (Control)');
%% Feeding
% % Select dataset
% % Control
% M1F = load('20190228_M1F_sortedspikes.mat');
% M1U = load('20190228_M1U_sortedspikes.mat');
% M1F = struct2cell(M1F); 
% M1U = struct2cell(M1U);
% Units = {[M1F; M1U]};
% Units = Units{1};
% % Units = M1U;

% S1F = load('20190228_S1F_sortedspikes.mat');
% S1U = load('20190228_S1U_sortedspikes.mat');
% S1F = struct2cell(S1F); 
% S1U = struct2cell(S1U);
% Units = {[S1F; S1U]};
% Units = Units{1};
% % Units = S1U;

% M1F = load('20190509_M1F_sortedspikes.mat');
% M1U = load('20190509_M1U_sortedspikes.mat');
% M1F = struct2cell(M1F); 
% M1U = struct2cell(M1U);
% Units = {[M1F; M1U]};
% Units = Units{1};

S1F = load('20190509_S1F_sortedspikes.mat');
S1U = load('20190509_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = {[S1F; S1U]};
Units = Units{1};

%cut trials
conT = ControlTimes + 0.05;
spikecount = [];
timespikes = {};
meanfiring = [];
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(ControlTimes)
        inx=[];
        inx=find(temp >= ControlTimes(t)-0.2 & temp <= Cendtimes(t)+0.2);
        if inx == 0
            spikecount(t,u) = [];
        end
        spikecount(t,u) = length(inx);
        Ts = []; Ts = temp(inx);
        timespikes{t,u} = Ts-conT(t);
        meanfiring(t,u) = spikecount(t,u)/0.400; % over 400 ms
    end
%     meanfiring = ~isempty(meanfiring);
end
timespks = timespikes(ConI,:);

% Ranges of degrees
% directions = [-30,-20,-10,0,10,20,30];
directions = 1:9;
firingrate = [];
meanfr = [];
FRbyUnit = {};
for d = 1:(length(directions)-1)
    dir = ConI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit{n}(d,:) = firingrate(:,n);
    end
end
M1FCon = meanfr;

%Preferred Directions
pref = zeros(length(Units),1);
for unit = 1:length(Units)
    [p, d] = max(M1FCon(:,unit));
    pref(unit) = d;
end