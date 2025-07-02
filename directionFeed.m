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
indices = {}; times = {};
maxgapes = {};
tonguetip_x = {}; tonguetip_y = {}; tonguetip_z = {};
yaw = []; pitch = [];
a3 = []; a2 = []; f = [];
tt = {};
v = {};
%Control
for cycle = 1:length(cyclengths{1,1})
    trial = cyclemat{1,1}.Trialname{cycle};
    tempdata = Kinematics(1).Kinematics.Cranium.points{strcmp(trial,Kinematics(1).Kinematics.TrialNames)};
    allcycles = Controlcyc;

    temptimes = Kinematics(1).Kinematics.NeuralIndex{strcmp(trial,Kinematics(1).Kinematics.TrialNames)}(:,3);

    % Positions max gape to max gape
    tonguetip_x{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),28);
    tonguetip_y{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),29);
    tonguetip_z{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),30);
    
    % Instantaneous direction every 100 ms
    npz={[0 0 1] [0 -1 0]}; %reference plane
    fin = cyclengths{1,1}(cycle);
    tp=1:20:fin; %timepoints
    tp2=20:20:fin;
    if cyclengths{1,1}(cycle) < 40
        continue
    end
    for i = 1:length(tp2)
        v1 = [tonguetip_x{cycle}(tp(i)) tonguetip_z{cycle}(tp(i)) tonguetip_y{cycle}(tp(i))];
        v2 = [tonguetip_x{cycle}(tp2(i)) tonguetip_z{cycle}(tp2(i)) tonguetip_y{cycle}(tp2(i))];
        a3(cycle,i) = vecangle360(v1,v2,npz{1});
        a2(cycle,i) = vecangle360(v1,v2,npz{2});
        f(cycle,i) = tonguetip_x{cycle}(tp2(i)) - tonguetip_x{cycle}(tp(i));
        
        % Yaw and Pitch
        pt = v2-v1;
        yaw(cycle,i) = atan2(pt(3), pt(1));
        adj = sqrt(pt(1)^2 + pt(3)^2);
        pitch(cycle,i) = atan2(pt(2), adj);
        
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
        indices{cycle,i} = allcycles(cycle,1)+20*(i-1);
        times{cycle,i} = temptimes(indices{cycle,i});
    end
end
TT = tt'; TT = TT(:);
pos1 = find(~cellfun(@isempty,TT));
TT = TT(pos1);
ControlTimes = times'; ControlTimes = ControlTimes(:); ControlTimes = ControlTimes(pos1);
ControlTimes = cell2mat(ControlTimes);
Cendtimes = ControlTimes + 0.100;
% ControlTimes = Controltimes(pos1) - 0.05; %% shift back 50 ms
% Cendtimes = Controltimes(pos1) + 0.05;

% Distribution of directions (8 directions w/ A-P)
a3=a3'; A3=a3(:); a2=a2';A2=a2(:); f=f';F=f(:);
A3 = A3(pos1); A2 = A2(pos1); F = F(pos1);

cat = cell(length(A3),1);
vecM = cell(length(A3),1);
x = []; y = []; z = [];
for ang = 1:length(A3)
    if A3(ang) > 0 && A2(ang) > 0 && F(ang) > 0
        cat{ang} = '+++';
        vecM{ang} = [1,1,1]; % x,y,z
    end
    if A3(ang) < 0 && A2(ang) > 0 && F(ang) > 0
        cat{ang} = '-++';
        vecM{ang} = [1,1,-1];
    end
    if A3(ang) < 0 && A2(ang) < 0 && F(ang) > 0
        cat{ang} = '--+';
        vecM{ang} = [1,-1,-1];
    end
    if A3(ang) > 0 && A2(ang) < 0 && F(ang) > 0
        cat{ang} = '+-+';
        vecM{ang} = [1,-1,1];
    end    

    if A3(ang) > 0 && A2(ang) > 0 && F(ang) < 0
        cat{ang} = '++-';
        vecM{ang} = [-1,1,1];
    end
    if A3(ang) < 0 && A2(ang) > 0 && F(ang) < 0
        cat{ang} = '-+-';
        vecM{ang} = [-1,1,-1];
    end
    if A3(ang) < 0 && A2(ang) < 0 && F(ang) < 0
        cat{ang} = '---';
        vecM{ang} = [-1,-1,-1];
    end
    if A3(ang) > 0 && A2(ang) < 0 && F(ang) < 0
        cat{ang} = '+--';
        vecM{ang} = [-1,-1,1];
    end

    % Direction angles
    pt = vecM{ang};
    x(ang) = atan2(sqrt(pt(2)^2+pt(3)^2),pt(1));
    y(ang) = atan2(sqrt(pt(1)^2+pt(3)^2),pt(2));
    z(ang) = atan2(sqrt(pt(1)^2+pt(2)^2),pt(3));
end
x = x'; y = y'; z = z';

category = {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'};
% category = {'-++', '+++', '+-+','---','+--'}; %% swallows
catind = {}; mag = {};
% figure;
% t = tiledlayout(4,2,'TileSpacing','Compact');
for c = 1:length(category)
    catind{c} = find(contains(cat,category{c}));
    mag{c} = abs(A3(catind{c})); %angle

%     nexttile
%     histogram(mag{c});
%     title(category{c});
%     ylim([0 100]); %400
end
% title(t,'Distribution of directions')
% xlabel(t,'Angle (°)')

x = x'; y = y'; z = z';

% figure;hist(A3)
% figure;hist(A2)

yaw = yaw'; pitch = pitch';
Y = yaw(:); Y = Y(pos1); Yaw = circ_rad2ang(Y);
P = pitch(:); P = P(pos1); Pitch = circ_rad2ang(P);
% figure;hist(Yaw)
% figure;hist(Pitch)

% cosine directions
Ang1 = x(:);
Ang2 = y(:);
Ang3 = z(:);
mx = cos(Ang1); my = cos(Ang2); mz = cos(Ang3);

% ControlDir = A3;%mag{1};
% ControlInd = catind{1};
ControlDir = Yaw;
% ControlDir = Xcomp;
% figure;hist(A2);

% Trajectory Mat
ControlTraj = TT; %[tonguetip_x; tonguetip_y; tonguetip_z];


% Nerve Block
maxgapes = {};
tonguetip_x = {}; tonguetip_y = {}; tonguetip_z = {};
indices = {}; times = {};
a3 = []; a2 = []; f = [];
yaw = []; pitch = [];
tt = {};
for cycle = 1:length(cyclengths{1,2})
    trial = cyclemat{1,2}.Trialname{cycle};
    tempdata = Kinematics(2).Kinematics.Cranium.points{strcmp(trial,Kinematics(2).Kinematics.TrialNames)};
    allcycles = NBcyc;

    temptimes = Kinematics(2).Kinematics.NeuralIndex{strcmp(trial,Kinematics(2).Kinematics.TrialNames)}(:,3);

    % Positions max gape to max gape
    tonguetip_x{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),28);
    tonguetip_y{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),29);
    tonguetip_z{cycle} = tempdata(allcycles(cycle,1):allcycles(cycle,3),30);
    
    % Instantaneous direction every 100 ms
    npz={[0 0 1] [0 -1 0]};  %reference plane
    fin = cyclengths{1,2}(cycle);
    tp=1:20:fin; %timepoints
    tp2=20:20:fin;
    if cyclengths{1,2}(cycle) < 40
        continue
    end
    for i = 1:length(tp2)
        v1 = [tonguetip_x{cycle}(tp(i)) tonguetip_z{cycle}(tp(i)) tonguetip_y{cycle}(tp(i))];
        v2 = [tonguetip_x{cycle}(tp2(i)) tonguetip_z{cycle}(tp2(i)) tonguetip_y{cycle}(tp2(i))];
        a3(cycle,i) = vecangle360(v1,v2,npz{1});
        a2(cycle,i) = vecangle360(v1,v2,npz{2});
        f(cycle,i) = tonguetip_x{cycle}(tp2(i)) - tonguetip_x{cycle}(tp(i));

        pt = v2-v1;
        yaw(cycle,i) = atan2(pt(3), pt(1));
        adj = sqrt(pt(1)^2 + pt(3)^2);
        pitch(cycle,i) = atan2(pt(2), adj);
        
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
        indices{cycle,i} = allcycles(cycle,1)+20*(i-1);
        times{cycle,i} = temptimes(indices{cycle,i});
    end
end
TT = tt';
TT = TT(:);
pos2 = find(~cellfun(@isempty,TT));
TT = TT(pos2);
NBTimes = times'; NBTimes = NBTimes(:); NBTimes = NBTimes(pos2);
NBTimes = cell2mat(NBTimes);
NBendtimes = NBTimes + 0.100;
% NBTimes = NBtimes(pos2) - 0.05; %% shift back 50 ms
% NBendtimes = NBtimes(pos2) + 0.05;

% Distribution of directions
% A=a(:);
% A = A(pos2);
% figure;hist(A)
% Y = yaw(:);
% Y = Y(pos2);
% P = pitch(:);
% P = P(pos2);
% Yaw = circ_rad2ang(Y);
% Pitch = circ_rad2ang(P);
% figure;hist(Yaw)
% figure;hist(Pitch)

% Distribution of directions
a3=a3'; A3=a3(:); a2=a2';A2=a2(:); f=f';F=f(:);
A3 = A3(pos2); A2 = A2(pos2); F = F(pos2);

cat2 = cell(length(A3),1);
vecM2 = cell(length(A3),1);
x = []; y = []; z = [];
for ang = 1:length(A3)
    if A3(ang) > 0 && A2(ang) > 0 && F(ang) > 0
        cat2{ang} = '+++';
        vecM2{ang} = [1,1,1]; % x,y,z
    end
    if A3(ang) < 0 && A2(ang) > 0 && F(ang) > 0
        cat2{ang} = '-++';
        vecM2{ang} = [1,1,-1];
    end
    if A3(ang) < 0 && A2(ang) < 0 && F(ang) > 0
        cat2{ang} = '--+';
        vecM2{ang} = [1,-1,-1];
    end
    if A3(ang) > 0 && A2(ang) < 0 && F(ang) > 0
        cat2{ang} = '+-+';
        vecM2{ang} = [1,-1,1];
    end
    

    if A3(ang) > 0 && A2(ang) > 0 && F(ang) < 0
        cat2{ang} = '++-';
        vecM2{ang} = [-1,1,1];
    end
    if A3(ang) < 0 && A2(ang) > 0 && F(ang) < 0
        cat2{ang} = '-+-';
        vecM2{ang} = [-1,1,-1];
    end
    if A3(ang) < 0 && A2(ang) < 0 && F(ang) < 0
        cat2{ang} = '---';
        vecM2{ang} = [-1,-1,-1];
    end
    if A3(ang) > 0 && A2(ang) < 0 && F(ang) < 0
        cat2{ang} = '+--';
        vecM2{ang} = [-1,-1,1];
    end

    % Direction angles
    pt = vecM2{ang};
    x(ang) = atan2(sqrt(pt(2)^2+pt(3)^2),pt(1));
    y(ang) = atan2(sqrt(pt(1)^2+pt(3)^2),pt(2));
    z(ang) = atan2(sqrt(pt(1)^2+pt(2)^2),pt(3));

end

category2 = {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'};
catind2 = {};
mag2 = {};
% figure;
% t = tiledlayout(4,2,'TileSpacing','Compact');
for c = 1:length(category2)
    catind2{c} = find(contains(cat2,category2{c}));
    mag2{c} = abs(A3(catind2{c})); %angle

%     nexttile
%     histogram(mag2{c});
%     title(category2{c});
%     ylim([0 100]); %150
end
% title(t,'Distribution of directions')
% xlabel(t,'Angle (°)')

Ang1 = x(:);
Ang2 = y(:);
Ang3 = z(:);
mx2 = cos(Ang1); my2 = cos(Ang2); mz2 = cos(Ang3);

NBDir = A3; %mag2{1}
% NBInd = catind2{1};

% NBDir = A;
% NBDir = Yaw;

% Trajectory Mat
NBTraj = TT; %[tonguetip_x; tonguetip_y; tonguetip_z];
%% Sample Trials
% Pick equal number of trials from each range of directions

% Control
% directions = -30:20:30;%0:5:30; 10
directions = [-40,-20, -10,10, 20,40];
% directions = [-35,-25, -5,5, 25,35];
% sample = [];
% sampleidx = [];
% n = 45; % number of trials for each direction (45)
% for d = 1:(length(directions)-1)
%     dir = ControlDir(ControlDir > directions(d) & ControlDir < directions(d+1));
%     diridx = find(ControlDir > directions(d) & ControlDir < directions(d+1)); %ControlInd before find for ang
%         m = randperm(length(dir),n);
%         sample(d,:) = dir(m(1:n)); 
%         sampleidx(d,:) = diridx(m(1:n));
% 
% %         sample{d} = dir; 
% %         sampleidx{d} = diridx;
% end
% % ConI = sampleidx([1,3,5],:);

% 3D 4/8 directions
sample = [];
sampleidx = [];
n = 80; %% number of trials per direction 80 %60 swallows
for d = 1:length(catind)
    dir = mag{d};
    diridx = catind{d};
        m = randperm(length(catind{d}),n);
        sample(d,:) = dir(m(1:n)); 
        sampleidx(d,:) = diridx(m(1:n));
end

Con = sample(:);
ConI = sampleidx(:,:);
% figure;hist(Con)

% Nerve Block
% directions = -30:10:30; %0:5:30;
% directions = [-40,-20, -10,10, 20,40];
% directions = [-35,-25, -5,5, 25,35];
% sample2 = [];
% sampleidx2 = [];
% n = 45; %45
% for d = 1:(length(directions)-1)
%     dir = NBDir(NBDir > directions(d) & NBDir < directions(d+1));
%     diridx = find(NBDir > directions(d) & NBDir < directions(d+1));
%         m = randperm(length(dir),n);
%         sample2(d,:) = dir(m(1:n));   
%         sampleidx2(d,:) = diridx(m(1:n));
% end
% NBI = sampleidx2([1,3,5],:);

% 3D 4/8 directions
sample2 = [];
sampleidx2 = [];
n = 80;
for d = 1:length(catind2)
    dir = mag2{d};
    diridx = catind2{d};
        m = randperm(length(catind2{d}),n);
        sample2(d,:) = dir(m(1:n)); 
        sampleidx2(d,:) = diridx(m(1:n));
end

NB = sample2(:);
NBI = sampleidx2(:,:);
% figure;hist(NB)

%% MIo
% Control
M1F = load('20190509_M1F_sortedspikes.mat'); %% paths to neural data
M1U = load('20190509_M1U_sortedspikes.mat');
M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
Units = [M1F; M1U];
% Units = M1F; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
Counts = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(ControlTimes)

        % 1ms bins for FA %
        edges = ControlTimes(t) : 0.001 : (Cendtimes(t) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        Counts{t}(u,:) = allCounts;

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
        timespikes{t,u} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
% directions = -30:20:30;%0:5:30;
directions = 1:9; %6 for swallows
% directions = 1:4;
firingrate = [];
meanfr = [];
FRbyUnit = {};
SpikesbyDir = {};
for d = 1:(length(directions)-1)
%     dir = ConI(d,:);
%     ncycles = length(dir);
%     cycles = meanfiring(dir,:);
    dirAll = catind{d}'; %all
    ncycles = length(dirAll); %all
    cycles = meanfiring(dirAll,:); %all
    SpikesbyDir{d} = spikecount(dirAll,:); %all
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
%         FRbyUnit{n}(d,:) = firingrate(:,n);
        FRbyUnit{n}{d} = firingrate(:,n); %all
    end
end
M1Con = meanfr;

% % Tuning Curves
% yFit6 = [];
% angle = 1:8; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = M1Con(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
%     yFit6(:,unit) = cos_fun(p,angle);
%     m(unit) = max(yFit6(:,unit));
% end

% % Plot mean firing rate tuning curves
% angle = 1:(length(directions)-1); % x-axis
% % angle = 1:12;
% for u = 12 %length(Units)
%     figure;
%     axes1 = axes;
%     hold(axes1,'on');
%     plot(angle,M1Con(:,u),'LineWidth',1.5,'Color','r') % actual data
%     hold on
% %     plot(angle,yFit6(:,u),'r-.') % cosine fit
% %     legend('Actual Data','Cosine Fit')
% end
% % ylabel('Average Firing Rate');
% % xlabel('3D angle (°)');
% title('MIo');
% xlim(axes1,[0.5 3.5]);
% box(axes1,'on');
% hold(axes1,'off');
% set(axes1,'FontSize',18,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0','0:5','5:10','10:15','15:20','20:25','25:30'});
%     {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'});

%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0','0:5','5:10','10:15','15:20','20:25','25:30'});
% {'0:5','5:10','10:15','15:20','20:25','25:30'});
%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0'});

% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0','0:5','5:10','10:15','15:20','20:25','25:30'});
% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32],'XTickLabel',...
%     {'-160:-150','-150:-140','-140:-130','-130:-120','-120:-110','-110:-100','-100:-90','-90:-80','-80:-70','-70:-60','-60:-50','-50:-40','-40:-30','-30:-20','-20:-10','-10:0','0:10','10:20','20:30','30:40','40:50','50:60','60:70','70:80','80:90','90:100','100:110','110:120','120:130','130:140','140:150','150:160'});


% Nerve Block
M1F = load('20190510_M1F_sortedspikes.mat'); %% paths to neural data
M1U = load('20190510_M1U_sortedspikes.mat');
M1F = struct2cell(M1F); 
M1U = struct2cell(M1U);
Units = [M1F; M1U];
% Units = M1F; %% run individual area

spikecount = []; timespikes = {}; meanfiring = [];
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(NBTimes)

        % 1ms bins for FA %
        edges = NBTimes(t) : 0.001 : (NBendtimes(t) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        NBCounts{t}(u,:) = allCounts;

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
        timespikes{t,u} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
% directions = -30:10:30;
% directions = 1:9; 
directions = 1:4;
firingrate = [];
meanfr = [];
FRbyUnit2 = {};
SpikesbyDir = {};
for d = 1:(length(directions)-1)
    dir = NBI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
%     dirAll = catind2{d}'; %all
%     ncycles = length(dirAll); %all
%     cycles = meanfiring(dirAll,:); %all
%     SpikesbyDir{d} = spikecount(dirAll,:); %all
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit2{n}(d,:) = firingrate(:,n);
%         FRbyUnit2{n}{d} = firingrate(:,n); %all
    end
end
M1NB = meanfr;

% % Tuning Curves
% yFit = [];
% yFit8 = [];
% angle = 1:8; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = M1NB(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
%     p = nlinfit(angle,y,cos_fun,[1 1 0]);
% %     yFit(:,unit) = cos_fun(p,[-30:30]);
%     yFit8(:,unit) = cos_fun(p,angle);
% %     m(unit) = max(yFit(:,unit));
% end

% % Plot mean firing rate tuning curves
% for u = 45 %length(Units)
%     figure;
%     axes1 = axes;
%     hold(axes1,'on');
%     plot(angle,M1NB(:,u)) % actual data
%     xlabel('Angle')
%     ylabel('Average Firing Rate')
%     hold on
% %     plot(angle,yFit6(:,u),'r-.') % cosine fit
%     legend('Actual Data','Cosine Fit')
% end
% ylabel('Average Firing Rate');
% xlabel('Angle (deg)');
% title('MIo');
% xlim(axes1,[0.5 4.5]);
% box(axes1,'on');
% hold(axes1,'off');
% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0','0:5','5:10','10:15','15:20','20:25','25:30'});

% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32],'XTickLabel',...
%     {'-160:-150','-150:-140','-140:-130','-130:-120','-120:-110','-110:-100','-100:-90','-90:-80','-80:-70','-70:-60','-60:-50','-50:-40','-40:-30','-30:-20','-20:-10','-10:0','0:10','10:20','20:30','30:40','40:50','50:60','60:70','70:80','80:90','90:100','100:110','110:120','120:130','130:140','140:150','150:160'});

%% SIo
% Control
S1F = load('20190509_S1F_sortedspikes.mat'); %% paths to neural data
S1U = load('20190509_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = [S1F; S1U];
% Units = S1F; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
Counts = {};
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(ControlTimes)

        % 1ms bins for FA %
        edges = ControlTimes(t) : 0.001 : (Cendtimes(t) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        Counts{t}(u,:) = allCounts;

        inx=[];
        inx=find(temp >= ControlTimes(t) & temp <= Cendtimes(t));
        spikecount(t, u) = length(inx);
        Ts = []; Ts = temp(inx);
        if ~isempty(Ts)
            Ts = Ts - temp(1);
        end
        timespikes{t,u} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
% directions = -30:20:30;
% directions = 1:9;%6;
firingrate = [];
meanfr = [];
FRbyUnit = {};
SpikesbyDir = {};
for d = 1:(length(directions)-1)
    dir = ConI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
%     dirAll = catind{d}';
%     ncycles = length(dirAll);
%     cycles = meanfiring(dirAll,:);
%     SpikesbyDir{d} = spikecount(dirAll,:);
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit{n}(d,:) = firingrate(:,n);
%         FRbyUnit{n}{d} = firingrate(:,n); %all
    end
end
S1Con = meanfr;

% % Tuning Curve
% yFit6 = [];
% angle = 1:32; % x-axis
% m = [];
% for unit = 1:length(Units)
%     y = S1Con(:,unit)';
%     cos_fun = @(p,theta) p(1)+p(2)*cos(theta/180-p(3));
% %     p = nlinfit(angle,y,cos_fun,[1 1 0]);
%     yFit6(:,unit) = cos_fun(p,angle);
% %     m(unit) = max(yFit6(:,unit));
% end

% % Plot mean firing rate tuning curves
% angle = 1:(length(directions)-1);
% for u = 2
%     figure;
%     axes1 = axes;
%     hold(axes1,'on');
%     plot(angle,S1Con(:,u),'LineWidth',1.5,'Color','r') % actual data
%     hold on
% %     plot([-30:30],yFit6(:,u),'r-.') % cosine fit
% %     legend('Actual Data','Cosine Fit')
% end
% % ylabel('Average Firing Rate');
% % xlabel('3D angle (°)');
% title('SIo');
% xlim(axes1,[0.5 4.5]);
% box(axes1,'on');
% hold(axes1,'off');
% set(axes1,'FontSize',18,'XTick',[1 2 3 4 5 6 7 8],'XTickLabel',...
%     {'-++', '+++', '--+', '+-+','-+-','++-','---','+--'});
% %         {'0:5','5:10','10:15','15:20','20:25','25:30'});

% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12],'XTickLabel',...
%     {'-30:-25','-25:-20','-20:-15','-15:-10','-10:-5','-5:0','0:5','5:10','10:15','15:20','20:25','25:30'});
% set(axes1,'FontSize',12,'XTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32],'XTickLabel',...
%     {'-160:-150','-150:-140','-140:-130','-130:-120','-120:-110','-110:-100','-100:-90','-90:-80','-80:-70','-70:-60','-60:-50','-50:-40','-40:-30','-30:-20','-20:-10','-10:0','0:10','10:20','20:30','30:40','40:50','50:60','60:70','70:80','80:90','90:100','100:110','110:120','120:130','130:140','140:150','150:160'});


% Nerve Block
S1F = load('20190510_S1F_sortedspikes.mat'); %% paths to neural data
S1U = load('20190510_S1U_sortedspikes.mat');
S1F = struct2cell(S1F); 
S1U = struct2cell(S1U);
Units = [S1F; S1U];
% Units = S1U; %% run individual area

spikecount = [];
timespikes = {};
meanfiring = [];
for u = 1:length(Units)
    temp = Units{u}.times;
    for t=1:length(NBTimes)

        % 1ms bins for FA %
        edges = NBTimes(t) : 0.001 : (NBendtimes(t) + 0.001); %
        allCounts = zeros(1, length(edges) - 1); %
        [counts, values] = histcounts(temp, edges); %
        allCounts = allCounts + counts; %
        NBCounts{t}(u,:) = allCounts;

        inx=[];
        inx=find(temp >= NBTimes(t) & temp <= NBendtimes(t));
        spikecount(t, u) = length(inx);
        Ts = []; Ts = temp(inx);
        if ~isempty(Ts)
            Ts = Ts - temp(1);
        end
        timespikes{t,u} = Ts;
        meanfiring(t,u) = spikecount(t,u)/0.100; %% 100 ms
    end
end

% Ranges of degrees
% directions = -30:10:30;
% directions = 1:9;
firingrate = [];
meanfr = [];
FRbyUnit2 = {};
SpikesbyDir = {};
for d = 1:(length(directions)-1)
    dir = NBI(d,:);
    ncycles = length(dir);
    cycles = meanfiring(dir,:);
%     dirAll = catind2{d}'; %all
%     ncycles = length(dirAll); %all
%     cycles = meanfiring(dirAll,:); %all
%     SpikesbyDir{d} = spikecount(dirAll,:); %all
    firingrate = cycles;
    for n = 1:length(Units)
        meanfr(d,n) = sum(firingrate(:,n))/ncycles;
        FRbyUnit2{n}(d,:) = firingrate(:,n);
%         FRbyUnit2{n}{d} = firingrate(:,n); %all
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
