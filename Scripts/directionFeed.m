%% Get directions and spike data (feeding)
%% Load Data
Control = load('20190509_Kinematics.mat');
NerveBlock = load('20190510_Kinematics.mat');
Kinematics = [Control NerveBlock];
% copy and paste paired dataset paths
%% Instantaneous 3D Directions
cyclemat = {Kinematics(1).Kinematics.GapeCycleInfo Kinematics(2).Kinematics.GapeCycleInfo};

% chews = find(contains(Kinematics(1).Kinematics.GapeCycleInfo.CycleType,'Chew'));
% cyclemat{1} = cyclemat{1}(chews,:); %% chews/swallows ONLY
% chews2 = find(contains(Kinematics(2).Kinematics.GapeCycleInfo.CycleType,'Chew'));
% cyclemat{2} = cyclemat{2}(chews2,:); %% chews/swallows ONLY

cyclengths = {cyclemat{1,1}.MaxGapeEnd - cyclemat{1,1}.MaxGapeStart cyclemat{1,2}.MaxGapeEnd - cyclemat{1,2}.MaxGapeStart};
Controlcyc = table2array(cyclemat{1,1}(:,9:11));
NBcyc = table2array(cyclemat{1,2}(:,9:11));
maxgapes = {};
tonguetip_x = {};
tonguetip_y = {};
tonguetip_z = {};
a = [];
tt = {};
%Control
for cycle = 1:length(cyclengths{1,1})
    trial = cyclemat{1,1}.Trialname{cycle};
    tempdata = Kinematics(1).Kinematics.Cranium.points{strcmp(trial,Kinematics(1).Kinematics.TrialNames)};
    allcycles = Controlcyc;

    % Positions max gape to max gape
    tonguetip_x{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),28);
    tonguetip_y{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),29);
    tonguetip_z{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),30);
    
    % Instantaneous direction every 100 ms
    npz=[0 0 1]; %reference plane
    fin = cyclengths{1,1}(cycle);
    tp=1:20:fin; %timepoints
    tp2=20:20:fin;
    if cyclengths{1,1}(cycle) < 40
        continue
    end
    for i = 1:length(tp2)
        v1 = [tonguetip_x{cycle}(tp(i)) tonguetip_z{cycle}(tp(i)) tonguetip_y{cycle}(tp(i))];
        v2 = [tonguetip_x{cycle}(tp2(i)) tonguetip_z{cycle}(tp2(i)) tonguetip_y{cycle}(tp2(i))];
        a(cycle,i) = vecangle360(v1,v2,npz);
        
%         % Visualization of trajectories
%         figure(i)
%         subplot(1,2,1)
%         plot3(tonguetip_x{cycle}(tp(i):tp2(i)),tonguetip_z{cycle}(tp(i):tp2(i)),tonguetip_y{cycle}(tp(i):tp2(i)),'k'); hold on;
%         plot3(tonguetip_x{cycle}(tp(i):tp2(i)),tonguetip_z{cycle}(tp(i):tp2(i)),tonguetip_y{cycle}(tp(i):tp2(i)),'r.')
%         plot3(tonguetip_x{cycle}(tp(i)),tonguetip_z{cycle}(tp(i)),tonguetip_y{cycle}(tp(i)),'bo')
%         plot3(tonguetip_x{cycle}(tp2(i)),tonguetip_z{cycle}(tp2(i)),tonguetip_y{cycle}(tp2(i)),'go')
%         axis square
%         
%         figure(i)
%         subplot(1,2,2)
%         plot(a(:,i),'b');hold on; %plot(a2,'r');plot(a3,'g');
%         axis square
        tt{cycle,i} = [tonguetip_x{cycle}(tp(i):tp2(i)) tonguetip_z{cycle}(tp(i):tp2(i)) tonguetip_y{cycle}(tp(i):tp2(i))];
   end
end
TT = tt(:);
pos1 = find(~cellfun(@isempty,TT));
TT = TT(pos1);

% Distribution of directions
A=a(:);
A = A(pos1);
% figure;hist(A)

ControlDir = A;

% Trajectory Mat
ControlTraj = TT; %[tonguetip_x; tonguetip_y; tonguetip_z];

% Nerve Block
maxgapes = {};
tonguetip_x = {};
tonguetip_y = {};
tonguetip_z = {};
a = [];
tt = {};
for cycle = 1:length(cyclengths{1,2})
    trial = cyclemat{1,2}.Trialname{cycle};
    tempdata = Kinematics(2).Kinematics.Cranium.points{strcmp(trial,Kinematics(2).Kinematics.TrialNames)};
    allcycles = NBcyc;

    % Positions max gape to max gape
    tonguetip_x{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),28);
    tonguetip_y{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),29);
    tonguetip_z{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),30);
    
    % Instantaneous direction every 100 ms
    npz=[0 0 1]; %reference plane
    fin = cyclengths{1,2}(cycle);
    tp=1:20:fin; %timepoints
    tp2=20:20:fin;
    if cyclengths{1,2}(cycle) < 40
        continue
    end
    for i = 1:length(tp2)
        v1 = [tonguetip_x{cycle}(tp(i)) tonguetip_z{cycle}(tp(i)) tonguetip_y{cycle}(tp(i))];
        v2 = [tonguetip_x{cycle}(tp2(i)) tonguetip_z{cycle}(tp2(i)) tonguetip_y{cycle}(tp2(i))];
        a(cycle,i) = vecangle360(v1,v2,npz);
        
%         % Visualization of trajectories
%         figure(i)
%         subplot(1,2,1)
%         plot3(tonguetip_x{cycle}(tp(i):tp2(i)),tonguetip_z{cycle}(tp(i):tp2(i)),tonguetip_y{cycle}(tp(i):tp2(i)),'k'); hold on;
%         plot3(tonguetip_x{cycle}(tp(i):tp2(i)),tonguetip_z{cycle}(tp(i):tp2(i)),tonguetip_y{cycle}(tp(i):tp2(i)),'r.')
%         plot3(tonguetip_x{cycle}(tp(i)),tonguetip_z{cycle}(tp(i)),tonguetip_y{cycle}(tp(i)),'bo')
%         plot3(tonguetip_x{cycle}(tp2(i)),tonguetip_z{cycle}(tp2(i)),tonguetip_y{cycle}(tp2(i)),'go')
%         axis square
%         
%         figure(i)
%         subplot(1,2,2)
%         plot(a(:,i),'b');hold on; %plot(a2,'r');plot(a3,'g');
%         axis square

        tt{cycle,i} = [tonguetip_x{cycle}(tp(i):tp2(i)) tonguetip_z{cycle}(tp(i):tp2(i)) tonguetip_y{cycle}(tp(i):tp2(i))];
   end
end
TT = tt(:);
pos2 = find(~cellfun(@isempty,TT));
TT = TT(pos2);

% Distribution of directions
A=a(:);
A = A(pos2);
% figure;hist(A)

NBDir = A;

% Trajectory Mat
NBTraj = TT; %[tonguetip_x; tonguetip_y; tonguetip_z];
%% Sample Trials
% Pick 45 trials from each range of directions
%Control
directions = [-30,-20,-10,0,10,20,30];
% directions = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]; %% for tuning curves
sample = [];
sampleidx = [];
all = {};
allidx = {};
for d = 1:6 %% 6 or 12
    dir = ControlDir(ControlDir > directions(d) & ControlDir < directions(d+1));
    diridx = find(ControlDir > directions(d) & ControlDir < directions(d+1));
        m = randperm(45);
        sample(d,:) = dir(m(1:45)); %% 45 or 20 for tuning curves
        sampleidx(d,:) = diridx(m(1:45));
end

Con = sample(:);
ConI = sampleidx(:,:);
% figure;hist(Con)

% Nerve Block
directions = [-30,-20,-10,0,10,20,30];
sample = [];
sampleidx = [];
for d = 1:6
    dir = NBDir(NBDir > directions(d) & NBDir < directions(d+1));
    diridx = find(NBDir > directions(d) & NBDir < directions(d+1));
        m = randperm(45);
        sample(d,:) = dir(m(1:45));   
        sampleidx(d,:) = diridx(m(1:45));
end

NB = sample(:);
NBI = sampleidx(:,:);
% figure;hist(NB)

%% Neural Times
% Get neural times of max gapes
% Control
maxtimes = {};
for cycle = 1:length(cyclengths{1,1})
    trial = cyclemat{1,1}.Trialname{cycle};
    tempdata = Control.Kinematics.NeuralIndex{strcmp(trial,Control.Kinematics.TrialNames)};
    tempdata = tempdata(:,3);
    size = length(1:20:cyclengths{1,1}(cycle));
    maxtimes{cycle,1} = tempdata(Controlcyc(cycle,1));
    for num = 1:size
        maxtimes{cycle,num+1} = maxtimes{cycle,num} + 3000;
        tf = cellfun('isempty',maxtimes);
        maxtimes(tf) = {0} ;
    end
end
ControlTimes = cell2mat(maxtimes(:));
ControlTimes = ControlTimes(pos1);
Cendtimes = ControlTimes + 0.100;
% ControlTimes = Controltimes(pos1) - 0.05; %% shift back 50 ms
% Cendtimes = Controltimes(pos1) + 0.05;

% Nerve Block
maxtimes = {};
for cycle = 1:length(cyclengths{1,2})
    trial = cyclemat{1,2}.Trialname{cycle};
    tempdata = NerveBlock.Kinematics.NeuralIndex{strcmp(trial,NerveBlock.Kinematics.TrialNames)};
    tempdata = tempdata(:,3);
    size = length(1:20:cyclengths{1,2}(cycle));
    maxtimes{cycle,1} = tempdata(NBcyc(cycle,1));
    for num = 1:size
        maxtimes{cycle,num+1} = maxtimes{cycle,num}+ 3000;
        tf = cellfun('isempty',maxtimes);
        maxtimes(tf) = {0} ;
    end
end
NBTimes = cell2mat(maxtimes(:));
NBTimes = NBTimes(pos2);
NBendtimes = NBTimes + 0.100;
% NBTimes = NBtimes(pos2) - 0.05; %% shift back 50 ms
% NBendtimes = NBtimes(pos2) + 0.05;
%% MIo
% Control
M1F = load('20190509_M1F_sortedspikes.mat'); %% paths to neural data
M1U = load('20190509_M1U_sortedspikes.mat');
M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
Units = {[M1F; M1U]};
Units = Units{1};
% Units = M1U; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(ControlTimes)
        inx=[];
        inx=find(temp >= ControlTimes(t) & temp <= Cendtimes(t));
        if inx == 0
            spikecount(t,u) = [];
        end
        spikecount(t,u) = length(inx);
        Ts = []; Ts = temp(inx);
        if ~isempty(Ts)
            Ts = Ts - temp(1);
        end
        timespikes{t} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
directions = [-30,-20,-10,0,10,20,30];
% directions = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30];
firingrate = [];
meanfr = [];
FRbyUnit = {};
for d = 1:6  %% number of directions 6 or 12
    dir = ConI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit{n}(d,:) = firingrate(:,n);
    end
end
M1Con = meanfr;

% % Tuning Curves
% yFit6 = [];
% angle = 1:12; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = M1Con(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
%     yFit6(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit6(:,unit));
% end
% 
% for u = 23 %length(Units)
%     figure;
%     axes1 = axes;
%     hold(axes1,'on');
%     plot(angle,M1Con(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
% %     plot(angle,yFit6(:,u),'r-.') % cosine fit
%     legend('Actual Data','Cosine Fit')
% end
% ylabel('Average Firing Rate');
% xlabel('Angle (deg)');
% title('MIo');
% xlim(axes1,[0.5 12.5]);
% box(axes1,'on');
% hold(axes1,'off');
% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0','0:5','5:10','10:15','15:20','20:25','25:30'});

% Nerve Block
M1F = load('20190510_M1F_sortedspikes.mat'); %% paths to neural data
M1U = load('20190510_M1U_sortedspikes.mat');
M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
Units = {[M1F; M1U]};
Units = Units{1};
% Units = M1U; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(NBTimes)
        inx=[];
        inx=find(temp >= NBTimes(t) & temp <= NBendtimes(t));
        if inx == 0
            spikecount(t,u) = [];
        end
        spikecount(t, u) = length(inx);
        Ts = []; Ts = temp(inx);
        if ~isempty(Ts)
            Ts = Ts - temp(1);
        end
        timespikes{t} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
directions = [-30,-20,-10,0,10,20,30];
firingrate = [];
meanfr = [];
FRbyUnit2 = {};
for d = 1:6
    dir = NBI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit2{n}(d,:) = firingrate(:,n);
    end
end
M1NB = meanfr;

% % Tuning Curves
% yFit = [];
% yFit6 = [];
% angle = [1:6]; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = M1NB(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
%     yFit(:,unit) = cos_fun(p,[-30:30]);
%     yFit8(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit(:,unit));
% end

% for u = 1:length(Units)
%     figure;
%     plot(angle,M1NB(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
%     %plot([-30:30],yFit(:,u),'r-.') % cosine fit
%     %legend('Actual Data','Cosine Fit')
% end
%% SIo
% Control
S1F = load('20190509_S1F_sortedspikes.mat'); %% paths to neural data
S1U = load('20190509_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = {[S1F; S1U]};
Units = Units{1};
% Units = S1U; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(ControlTimes)
        inx=[];
        inx=find(temp >= ControlTimes(t) & temp <= Cendtimes(t));
        spikecount(t, u) = length(inx);
        Ts = []; Ts = temp(inx);
        if ~isempty(Ts)
            Ts = Ts - temp(1);
        end
        timespikes{t} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
directions = [-30,-20,-10,0,10,20,30];
% directions = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30];
firingrate = [];
meanfr = [];
FRbyUnit = {};
for d = 1:6  %% number of directions 6 or 12
    dir = ConI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit{n}(d,:) = firingrate(:,n);
    end
end
S1Con = meanfr;

% % Tuning Curve
% yFit6 = [];
% angle = 1:12; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = S1Con(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
%     yFit6(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit6(:,unit));
% end
% 
% for u = 2
%     figure;
%     axes1 = axes;
%     hold(axes1,'on');
%     plot(angle,S1Con(:,u),'LineWidth',1.5,'Color','r') % actual data
%     xlabel('Angle (deg)')
%     ylabel('Average Firing Rate')
%     hold on
% %     plot([-30:30],yFit6(:,u),'r-.') % cosine fit
% %     legend('Actual Data','Cosine Fit')
% end
% ylabel('Average Firing Rate');
% xlabel('Angle (deg)');
% title('SIo');
% xlim(axes1,[0.5 12.5]);
% box(axes1,'on');
% hold(axes1,'off');
% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0','0:5','5:10','10:15','15:20','20:25','25:30'});

% Nerve Block
S1F = load('20190510_S1F_sortedspikes.mat'); %% paths to neural data
S1U = load('20190510_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = {[S1F; S1U]};
Units = Units{1};
% Units = S1U; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(NBTimes)
        inx=[];
        inx=find(temp >= NBTimes(t) & temp <= NBendtimes(t));
        spikecount(t, u) = length(inx);
        Ts = []; Ts = temp(inx);
        if ~isempty(Ts)
            Ts = Ts - temp(1);
        end
        timespikes{t} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
directions = [-30,-20,-10,0,10,20,30];
firingrate = [];
meanfr = [];
FRbyUnit2 = {};
for d = 1:6
    dir = NBI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit2{n}(d,:) = firingrate(:,n);
    end
end
S1NB = meanfr;

% % Tuning Curves
% yFit6 = [];
% angle = [1:6]; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = S1NB(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
%     yFit6(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit(:,unit));
% end

% for u = 1:length(Units)
%     figure;
%     plot(angle,S1NB(:,u)) % actual data
%     xlabel('Angle (deg)')
%     ylabel('Average Firing Rate')
%     hold on
%     %plot(angle,yFit6(:,u),'r-.') % cosine fit
%    legend('Actual Data','Cosine Fit')
% end