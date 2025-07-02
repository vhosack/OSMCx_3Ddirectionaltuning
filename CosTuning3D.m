%% 3D cosine tuning
% components for each direction
% multiple linear regression, check fit
% calculate and plot 3d preferred directions
% directional index + statistical tests
%% Get components
% Control
Mx = mx(sampleidx);
Mx = Mx(:);
My = my(sampleidx);
My = My(:);
Mz = mz(sampleidx);
Mz = Mz(:);

% % Nerve Block
% Mx = mx2(sampleidx2);
% Mx = Mx(:);
% My = my2(sampleidx2);
% My = My(:);
% Mz = mz2(sampleidx2);
% Mz = Mz(:);
%% Create models and check for fit
tuningmdl = {};
tuningtbl = {};
p = [];
for unit = 1:height(Tuned)   % Tuned2 if doing
    firing = Tuned{unit}(:); %   nerve block
    tuningtbl{unit} = table(Mx,My,Mz,firing);
    tuningmdl{unit} = fitlm(tuningtbl{unit});
    p(unit) = coefTest(tuningmdl{unit});
end

fits = p < 0.05;
percentFits = mean(fits)*100;

%% Get preferred directions
fitmdls = tuningmdl(fits); % pull only neurons that 
fittbls = tuningtbl(fits); %   fit the tuning function

figure;
Id = [];
pdir = {};
cx = []; cy = []; cz = [];
pds = [];
for fm = 1:length(fitmdls)
    b = tuningmdl{fm}.Coefficients.Estimate(1);
    bx = tuningmdl{fm}.Coefficients.Estimate(2);
    by = tuningmdl{fm}.Coefficients.Estimate(3);
    bz = tuningmdl{fm}.Coefficients.Estimate(4);
    k = sqrt((bx^2)+(by^2)+(bz^2));
    Id(fm) = k/b;
    cx(fm) = bx/k; cy(fm) = by/k; cz(fm) = bz/k;
    pdir{fm} = [cx(fm),cy(fm),cz(fm)];
    xp = acos(cx(fm)); yp = acos(cy(fm)); zp = acos(cz(fm));

    % colors for M1
    c_map = parula;        % pick your map
    pds(fm) = pdir{fm}(1);

%     if cx < 0
        scatter3(pdir{fm}(1),pdir{fm}(3),pdir{fm}(2),30,pdir{fm}(1),'sq','filled');
% quiver3(0,0,0,pdir{fm}(1),pdir{fm}(3),pdir{fm}(2));%,pdir{fm}(1));
%     elseif cx > 0
%         scatter3(pdir{fm}(1),pdir{fm}(3),pdir{fm}(2),'+');
%     end
    hold on;
end
% colormap autumn;
colorbar;
grid off;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); clim([-1, 1]);
set(gca,'FontSize',12);
title('3D Preferred Directions', 'FontSize',14);
xlabel('Posterior-Anterior');
ylabel('Left-Right');
zlabel('Inferior-Superior');
xline(0, '--k', 'LineWidth', 1);
yline(0, '--k', 'LineWidth', 1);
plot3([0 0], [0 0], zlim, '--k', 'LineWidth', 1);
% view([90 0]);
% plotCircle3D([0 0 0],[1 0 0],0.33);
% plotCircle3D([0 0 0],[1 0 0],0.66);
% plotCircle3D([0 0 0],[1 0 0],1);
% s = sphere(10);

%% plot second dataset and separate color bar (optional)
c_dat = pds;         % pick the data you want to use to define color
pnt_colors = interp1(linspace(min(c_dat), max(c_dat), size(c_map, 1)),...
    c_map, c_dat(:));
for fm = 1:length(fitmdls)
    scatter3(pdir{fm}(1),pdir{fm}(3),pdir{fm}(2),40,pnt_colors(fm,:),'x');
end
colorbar

%% plot costheta
figure;
dirgroups = 1:90:721;
fm = 23;
dM = [];
costheta = [];
for dir = 1:8
    dM(dir) = mean(fittbls{fm}.firing(dirgroups(dir):dirgroups(dir+1)-1));
    costheta(dir) = ((dM(dir) - b)/k);
    scatter(costheta(dir),dM(dir),20,'o','filled','k');
    hold on;
end
plot(costheta,dM,'k','LineWidth',1);

%% plot Id

figure;
% histogram(Id,'Normalization','probability');
histogram(Id,'FaceAlpha',0.8,'FaceColor',[0.494117647058824 0.184313725490196 0.556862745098039],...
    'Normalization','probability','BinMethod','auto');

ylabel('Proportion of neurons');
xlabel('Directional Index');
title('SIo');

set(gca,'FontSize',14,'YGrid','on');

%% Mean PD + variance
R = sqrt(sum(cx)^2 + sum(cy)^2 + sum(cz)^2)
meanPD = [sum(cx)/R, sum(cy)/R, sum(cz)/R]
sphvar = (length(fitmdls)-R)/length(fitmdls)

%% Rayleigh test
ax = acos(cx); ay = acos(cy); az = acos(cz);
[p1,z1] = circ_rtest(ax)
[p2,z2] = circ_rtest(ay)
[p3,z3] = circ_rtest(az)

%% Concentration test
% axM = ax; ayM = ay; azM = az;
axS = ax; ayS = ay; azS = az;
%%
[pval, f] = circ_ktest(axM,axS)
[pval, f] = circ_ktest(ayM,ayS)
[pval, f] = circ_ktest(azM,azS)
