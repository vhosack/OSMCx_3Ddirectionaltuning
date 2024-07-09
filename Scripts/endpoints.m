%% Get endpoints, distance from mean, failed cycles
%% Load Data
Before = load('20181210_Kinematics_Before.mat');
After = load('20181210_Kinematics_After.mat');
%% Get trajectories (Before)
right = find(Before.Kinematics.SpoutNumber == 2);
middle = find(Before.Kinematics.SpoutNumber == 1);
left = find(Before.Kinematics.SpoutNumber == 3);

Right = Before.Kinematics.Cranium.points(right);
Middle = Before.Kinematics.Cranium.points(middle);
Left = Before.Kinematics.Cranium.points(left);

ttptcols = find(contains(Before.Kinematics.ColumnNames.points,'AnteriorM_'));

% Right trials
ntrials = length(Right);
ttpt_r = {};
cycles_r = {};
Cycles_r = {};
endpoints_r = {};
begpoints_r = {};
for i = 1:ntrials
    ttpt_r{i} = Right{i}(:,ttptcols);
    npeaks = length(ttpt_r{i});
    stretch = ttpt_r{i}(:,2);
    [~,minstretch_r] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_r] = findpeaks(stretch,'MinPeakProminence',10);
    Mins = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    mins = Mins(minstretch_r,:);
    Maxs = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    maxs = Maxs(maxstretch_r,:);
    minlength = min(length(minstretch_r), length(maxstretch_r));
    minstretch_r = minstretch_r(1:minlength);
    maxstretch_r = maxstretch_r(1:minlength);
    stretchlength_r{i} = maxstretch_r - minstretch_r;
    for k = 1:minlength
        cycles_r{i,k}(:,1) = ttpt_r{i}(minstretch_r(k):maxstretch_r(k),1);
        cycles_r{i,k}(:,2) = ttpt_r{i}(minstretch_r(k):maxstretch_r(k),2);
        cycles_r{i,k}(:,3) = ttpt_r{i}(minstretch_r(k):maxstretch_r(k),3);
    end
    Cycles_r{i} = cycles_r(~cellfun(@isempty,cycles_r))';

    % Get endpoints
    cyc = length(Cycles_r{i});
    for l = 1:cyc
        begs = length(Cycles_r{i}{l});
        begpoints_r{i}(l,:) = Cycles_r{i}{l}(begs,:);
        endpoints_r{i}(l,:) = Cycles_r{i}{l}(1,:);
    end
end

% Middle trials
ntrials = length(Middle);
ttpt_m = {};
cycles_m = {};
Cycles_m = {};
endpoints_m = {};
begpoints_m = {};
for i = 1:ntrials
    ttpt_m{i} = Middle{i}(:,ttptcols);
    npeaks = length(ttpt_m{i});
    stretch = ttpt_m{i}(:,2);
    [~,minstretch_m] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_m] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_m);
    Mins = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    maxs = Maxs(maxstretch_m,:);
    minlength = min(length(minstretch_m), length(maxstretch_m));
    minstretch_m = minstretch_m(1:minlength);
    maxstretch_m = maxstretch_m(1:minlength);
    stretchlength_m{i} = maxstretch_m - minstretch_m;
    for k = 1:minlength
        cycles_m{i,k}(:,1) = ttpt_m{i}(minstretch_m(k):maxstretch_m(k),1);
        cycles_m{i,k}(:,2) = ttpt_m{i}(minstretch_m(k):maxstretch_m(k),2);
        cycles_m{i,k}(:,3) = ttpt_m{i}(minstretch_m(k):maxstretch_m(k),3);
    end
    Cycles_m{i} = cycles_m(~cellfun(@isempty,cycles_m))';

    % Get endpoints
    cyc = length(Cycles_m{i});
    for l = 1:cyc
        begs = length(Cycles_m{i}{l});
        begpoints_m{i}(l,:) = Cycles_m{i}{l}(begs,:);
        endpoints_m{i}(l,:) = Cycles_m{i}{l}(1,:);
    end
end

% Left trials
ntrials = length(Left);
ttpt_l = {};
cycles_l = {};
Cycles_l = {};
endpoints_l = {};
begpoints_l = {};
for i = 1:ntrials
    ttpt_l{i} = Left{i}(:,ttptcols);
    npeaks = length(ttpt_l{i});
    stretch = ttpt_l{i}(:,2);
    [~,minstretch_l] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_l] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_l);
    Mins = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    maxs = Maxs(maxstretch_l,:);
    minlength = min(length(minstretch_l), length(maxstretch_l));
    minstretch_l = minstretch_l(1:minlength);
    maxstretch_l = maxstretch_l(1:minlength);
    stretchlength_l{i} = maxstretch_l - minstretch_l;
    for k = 1:minlength
        cycles_l{i,k}(:,1) = ttpt_l{i}(minstretch_l(k):maxstretch_l(k),1);
        cycles_l{i,k}(:,2) = ttpt_l{i}(minstretch_l(k):maxstretch_l(k),2);
        cycles_l{i,k}(:,3) = ttpt_l{i}(minstretch_l(k):maxstretch_l(k),3);
    end
    Cycles_l{i} = cycles_l(~cellfun(@isempty,cycles_l))';

    % Get endpoints
    cyc = length(Cycles_l{i});
    for l = 1:cyc
        begs = length(Cycles_l{i}{l});
        begpoints_l{i}(l,:) = Cycles_l{i}{l}(begs,:);
        endpoints_l{i}(l,:) = Cycles_l{i}{l}(1,:);
    end
end

begsbefore = {begpoints_l;begpoints_m;begpoints_r};
endsbefore = {endpoints_l;endpoints_m;endpoints_r};

%% Variance of endpoints (Control)
l = endsbefore{1};
lefts = cat(1,l{:});
lefts = rmmissing(lefts,1);
m = endsbefore{2};
middles = cat(1,m{:});
middles = rmmissing(middles,1);
r = endsbefore{3};
rights = cat(1,r{:});
rights = rmmissing(rights,1);
% var(lefts)

RyCon = {lefts middles rights};

% % Plot
% figure;
% plot3(lefts(:,3),lefts(:,2),lefts(:,1),'go');
% title('Control Endpoints')
% xlabel('Z')
% ylabel('Y')
% zlabel('X')
% hold on;
% plot3(middles(:,3),middles(:,2),middles(:,1),'bo');
% plot3(rights(:,3),rights(:,2),rights(:,1),'yo');
% zlim([40 65])
% xlim([-30 30])

%% Failed cycles
% Find mean endpoint
Lavg = mean(lefts,"omitnan");
Mavg = mean(middles,"omitnan");
Ravg = mean(rights,"omitnan");

% Calculate distance of each endpoint from mean
Ldist = [];
for i = 1:length(lefts)
    Ldist(i) = norm(lefts(i,:)-Lavg);
end
Lsd = std(Ldist,"omitnan");

Mdist = [];
for i = 1:length(middles)
    Mdist(i) = norm(middles(i,:)-Mavg);
end
Msd = std(Mdist,"omitnan");

Rdist = [];
for i = 1:length(rights)
    Rdist(i) = norm(rights(i,:)-Ravg);
end
Rsd = std(Rdist,"omitnan");

% If distance is > 2 std, mark as missed
Lmiss = sum(Ldist > 2*Lsd)/length(lefts);
Mmiss = sum(Mdist > 2*Msd)/length(middles);
Rmiss = sum(Rdist > 2*Rsd)/length(rights);

% Total proportion of missed cycles
AllMissed = (Lmiss+Mmiss+Rmiss)/3

% % Plot means
% plot3(Lavg(3),Lavg(2),Lavg(1),'ko',MarkerSize=12,MarkerFaceColor='k')
% plot3(Mavg(3),Mavg(2),Mavg(1),'ko',MarkerSize=12,MarkerFaceColor='k')
% plot3(Ravg(3),Ravg(2),Ravg(1),'ko',MarkerSize=12,MarkerFaceColor='k')

%% Get trajectories (After)
right = find(After.Kinematics.SpoutNumber == 2);
middle = find(After.Kinematics.SpoutNumber == 1);
left = find(After.Kinematics.SpoutNumber == 3);

Right = After.Kinematics.Cranium.points(right);
Middle = After.Kinematics.Cranium.points(middle);
Left = After.Kinematics.Cranium.points(left);

ttptcols = find(contains(After.Kinematics.ColumnNames.points,'AnteriorM_'));

% Right trials
ntrials = length(Right);
ttpt_r = {};
cycles_r = {};
Cycles_r = {};
endpoints_r = {};
begpoints_r = {};
for i = 1:ntrials
    ttpt_r{i} = Right{i}(:,ttptcols);
    npeaks = length(ttpt_r{i});
    stretch = ttpt_r{i}(:,2);
    [~,minstretch_r] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_r] = findpeaks(stretch,'MinPeakProminence',10);
    Mins = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    mins = Mins(minstretch_r,:);
    Maxs = [ttpt_r{i}(:,1), ttpt_r{i}(:,2), ttpt_r{i}(:,3)];
    maxs = Maxs(maxstretch_r,:);
    minlength = min(length(minstretch_r), length(maxstretch_r));
    minstretch_r = minstretch_r(1:minlength);
    maxstretch_r = maxstretch_r(1:minlength);
    stretchlength_r{i} = maxstretch_r - minstretch_r;
    for k = 1:minlength
        cycles_r{i,k}(:,1) = ttpt_r{i}(minstretch_r(k):maxstretch_r(k),1);
        cycles_r{i,k}(:,2) = ttpt_r{i}(minstretch_r(k):maxstretch_r(k),2);
        cycles_r{i,k}(:,3) = ttpt_r{i}(minstretch_r(k):maxstretch_r(k),3);
    end
    Cycles_r{i} = cycles_r(~cellfun(@isempty,cycles_r))';

    % Get endpoints
    cyc = length(Cycles_r{i});
    for l = 1:cyc
        begs = length(Cycles_r{i}{l});
        begpoints_r{i}(l,:) = Cycles_r{i}{l}(begs,:);
        endpoints_r{i}(l,:) = Cycles_r{i}{l}(1,:);
    end
end

% Middle trials
ntrials = length(Middle);
ttpt_m = {};
cycles_m = {};
Cycles_m = {};
endpoints_m = {};
begpoints_m = {};
times_m = {};
for i = 1:ntrials
    ttpt_m{i} = Middle{i}(:,ttptcols);
    npeaks = length(ttpt_m{i});
    stretch = ttpt_m{i}(:,2);
    [~,minstretch_m] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_m] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_m);
    Mins = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_m{i}(:,1), ttpt_m{i}(:,2), ttpt_m{i}(:,3)];
    maxs = Maxs(maxstretch_m,:);
    minlength = min(length(minstretch_m), length(maxstretch_m));
    minstretch_m = minstretch_m(1:minlength);
    maxstretch_m = maxstretch_m(1:minlength);
    stretchlength_m{i} = maxstretch_m - minstretch_m;
    for k = 1:minlength
        cycles_m{i,k}(:,1) = ttpt_m{i}(minstretch_m(k):maxstretch_m(k),1);
        cycles_m{i,k}(:,2) = ttpt_m{i}(minstretch_m(k):maxstretch_m(k),2);
        cycles_m{i,k}(:,3) = ttpt_m{i}(minstretch_m(k):maxstretch_m(k),3);
    end
    Cycles_m{i} = cycles_m(~cellfun(@isempty,cycles_m))';

    % Get endpoints
    cyc = length(Cycles_m{i});
    for l = 1:cyc
        begs = length(Cycles_m{i}{l});
        begpoints_m{i}(l,:) = Cycles_m{i}{l}(begs,:);
        endpoints_m{i}(l,:) = Cycles_m{i}{l}(1,:);
    end
end

% Left trials
ntrials = length(Left);
ttpt_l = {};
cycles_l = {};
Cycles_l = {};
endpoints_l = {};
begpoints_l = {};
for i = 1:ntrials
    ttpt_l{i} = Left{i}(:,ttptcols);
    npeaks = length(ttpt_l{i});
    stretch = ttpt_l{i}(:,2);
    [~,minstretch_l] = findpeaks(stretch*-1,'MinPeakProminence',10);
    [~,maxstretch_l] = findpeaks(stretch,'MinPeakProminence',10);
    npeaks = length(minstretch_l);
    Mins = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    mins = Mins(minstretch_m,:);
    Maxs = [ttpt_l{i}(:,1), ttpt_l{i}(:,2), ttpt_l{i}(:,3)];
    maxs = Maxs(maxstretch_l,:);
    minlength = min(length(minstretch_l), length(maxstretch_l));
    minstretch_l = minstretch_l(1:minlength);
    maxstretch_l = maxstretch_l(1:minlength);
    stretchlength_l{i} = maxstretch_l - minstretch_l;
    for k = 1:minlength
        cycles_l{i,k}(:,1) = ttpt_l{i}(minstretch_l(k):maxstretch_l(k),1);
        cycles_l{i,k}(:,2) = ttpt_l{i}(minstretch_l(k):maxstretch_l(k),2);
        cycles_l{i,k}(:,3) = ttpt_l{i}(minstretch_l(k):maxstretch_l(k),3);
    end
    Cycles_l{i} = cycles_l(~cellfun(@isempty,cycles_l))';

    % Get endpoints
    cyc = length(Cycles_l{i});
    for l = 1:cyc
        begs = length(Cycles_l{i}{l});
        begpoints_l{i}(l,:) = Cycles_l{i}{l}(begs,:);
        endpoints_l{i}(l,:) = Cycles_l{i}{l}(1,:);
    end
end

endsafter = {endpoints_l;endpoints_m;endpoints_r};
begsafter = {begpoints_l;begpoints_m;begpoints_r};
%% Variance of endpoints (Nerve Block)
l = endsafter{1};
lefts = cat(1,l{:});
lefts = rmmissing(lefts,1);
m = endsafter{2};
middles = cat(1,m{:});
middles = rmmissing(middles,1);
r = endsafter{3};
rights = cat(1,r{:});
rights = rmmissing(rights,1);
% var(lefts)

RyNB = {lefts middles rights};

% % Plot
% figure;
% plot3(lefts(:,3),lefts(:,2),lefts(:,1),'go');
% title('Nerve Block Endpoints')
% xlabel('Z')
% ylabel('Y')
% zlabel('X')
% hold on;
% plot3(middles(:,3),middles(:,2),middles(:,1),'bo');
% plot3(rights(:,3),rights(:,2),rights(:,1),'yo');
% zlim([40 65])
% xlim([-30 30])

%% Failed cycles
% Find mean endpoint
Lavg = mean(lefts,"omitnan");
Mavg = mean(middles,"omitnan");
Ravg = mean(rights,"omitnan");

% Calculate distance of each endpoint from mean
Ldist = [];
for i = 1:length(lefts)
    Ldist(i) = norm(lefts(i,:)-Lavg);
end
Lsd = std(Ldist,"omitnan");

Mdist = [];
for i = 1:length(middles)
    Mdist(i) = norm(middles(i,:)-Mavg);
end
Msd = std(Mdist,"omitnan");

Rdist = [];
for i = 1:length(rights)
    Rdist(i) = norm(rights(i,:)-Ravg);
end
Rsd = std(Rdist,"omitnan");

% If distance is > 2 std, mark as missed
Lmiss = sum(Ldist > 2*Lsd)/length(lefts);
Mmiss = sum(Mdist > 2*Msd)/length(middles);
Rmiss = sum(Rdist > 2*Rsd)/length(rights);

% Total proportion of missed cycles
AllMissed = (Lmiss+Mmiss+Rmiss)/3

% % Plot means
% plot3(Lavg(3),Lavg(2),Lavg(1),'ko',MarkerSize=12,MarkerFaceColor='k')
% plot3(Mavg(3),Mavg(2),Mavg(1),'ko',MarkerSize=12,MarkerFaceColor='k')
% plot3(Ravg(3),Ravg(2),Ravg(1),'ko',MarkerSize=12,MarkerFaceColor='k')
