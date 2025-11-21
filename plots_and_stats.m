
% Sections à run séparément

%% Kinectomes

p=4;
load Kinect.mat
figure
tl=tiledlayout(2,3);
markers={'HEAD' 'C7' 'T10' 'CLAV' 'STRN' 'RBAK' 'LSHO' 'LELB' 'LFIN' 'RSHO' 'RELB' 'RFIN' ...
'LASI' 'RASI' 'LPSI' 'RPSI' 'LTHI' 'LKNE' 'LTIB' 'LANK' 'LHEE' 'LTOE' 'RTHI' 'RKNE' 'RTIB' 'RANK' 'RHEE' 'RTOE'};
cond=["Plat" "Medium" "High"];
for d=1:2
    for c=1:3
        nexttile
        h(c+(d-1)*3)=heatmap(markers,markers,Kinect{c,p}{end-1,d},'Colormap',parula,'ColorbarVisible','off'); % ,'GridVisible','off'
        if c==1 && d==1
            ylabel("Antéro-postérieure")
        elseif c==1 && d==2
            ylabel("Médio-latérale")
            xlabel(cond(1,c))
        elseif d==2
            xlabel(cond(1,c))
        end
    end
end

colorLims = vertcat(h.ColorLimits);
globalColorLim=[min(colorLims(:,1)),max(colorLims(:,2))];
set(h,'ColorLimits',globalColorLim)

ax=axes(tl,'visible','off','Colormap',h(1).Colormap,'CLim',globalColorLim);
cb=colorbar(ax);
cb.Layout.Tile = 'East';
ylabel(tl, "Direction",'FontWeight','bold')
xlabel(tl, "Condition",'FontWeight','bold')
tl.Padding = 'compact'; tl.TileSpacing = 'compact';

%% Plot MARP ou DV trois cond

paire=1; % 1 : Knee/Hip, 2 : Ankle/Knee
idms=1;  % 1 : MARP, 0 : DP

load PA_CRP.mat
load participants.mat
nbp=size(PA_CRP,2);

cond=["Even" "Medium" "High"];
group=["Adultes" "Adolescents" "Enfants" "Jeunes Enfants"];
ind=ones(3,4);
CRPmp=cell(3,4);
for p=1:nbp
    if isempty(PA_CRP{1,p})
        continue
    end
    idg=participants{p,3};
    for c=1:3
        CRPtemp(1,:)=PA_CRP{c,p}{end-idms,4}(paire,:);         % jambe gauche
        CRPtemp(2,:)=PA_CRP{c+3,p}{end-idms,4}(paire,:);       % jambe droite
        CRPmp{c,idg}(ind(c,idg),:)=mean(CRPtemp,1);
        ind(c,idg)=ind(c,idg)+1;
    end
end
CRPm=cell(3,4);
for i=1:3
    for j=1:4
        CRPm{i,j}=mean(CRPmp{i,j},1);
        CRPm{i,j}(2,:)=std(CRPmp{i,j},0,1);
    end
end

cms=colormap(turbo(4));

tl=tiledlayout(1,3);
for i=1:3
    nexttile
    for j=1:4
        hold on
        plot(CRPm{i,j}(1,:),'Color',cms(j,:),'LineWidth',2.5)
        f=fill([1:1:100 100:-1:1],[(CRPm{i,j}(1,:)+CRPm{i,j}(2,:)) fliplr((CRPm{i,j}(1,:)-CRPm{i,j}(2,:)))],'c');
        f.FaceColor=cms(j,:);
        f.EdgeColor='none';
        f.FaceAlpha=0.2;
    end
    title(cond(1,i),'FontWeight','bold')
    xlabel("% cycle")
    if i==1
        ylabel("MARP (°)")
    end
end
legend('Adultes','','Adolescents','','Enfants','','Jeunes Enfants','')
if paire==1
    title(tl,"Mean absolute relative phase - Knee/Hip")
else
    title(tl,"Mean absolute relative phase - Ankle/Knee")
end
tl.Padding = 'compact'; tl.TileSpacing = 'compact';

cms=colormap(parula(3));
figure
tl=tiledlayout(2,2);
for j=1:4
    nexttile
    for i=1:3
        hold on
        plot(CRPm{i,j}(1,:),'Color',cms(i,:),'LineWidth',2.5)
        f=fill([1:1:100 100:-1:1],[(CRPm{i,j}(1,:)+CRPm{i,j}(2,:)) fliplr((CRPm{i,j}(1,:)-CRPm{i,j}(2,:)))],'c');
        f.FaceColor=cms(i,:);
        f.EdgeColor='none';
        f.FaceAlpha=0.2;
    end
    title(group(1,j),'FontWeight','bold')
    if j==2
        legend('Even','','Medium','','High','')
    end
    ylabel("MARP (°)")
    xlabel("% cycle")
end
if paire==1
    title(tl,"Mean absolute relative phase - Knee/Hip")
else
    title(tl,"Mean absolute relative phase - Ankle/Knee")
end
tl.Padding = 'compact'; tl.TileSpacing = 'compact';

%% Post Hoc

paire=2; % 1 : Knee/Hip, 2 : Ankle/Knee
idms=1;  % 1 : MARP, 0 : DP

load PA_CRP.mat
load participants.mat
nbp=size(PA_CRP,2);

%--------------------------------------------------------------------------
% 
ind=1;
idp=0;
for p=1:nbp
    if isempty(PA_CRP{1,p})
        continue
    end
    idp=idp+1;
    for c=1:3
        crptmp=zeros(2,100);
        crptmp(1,:)=PA_CRP{c,p}{end-1,4}(paire,:);
        crptmp(2,:)=PA_CRP{c+3,p}{end-1,4}(paire,:);
        Y(ind,:)=mean(crptmp,1);
        A(ind,:)=participants{p,3};
        B(ind,:)=c;
        S(ind,:)=idp;
        ind=ind+1;
    end
end
% Pas de correction de sphéricité (equalvar dispo uniquement sur Python)
spm=spm1d.stats.normality.anova2onerm(Y,A,B,S);                             % Normality tests on the residuals
spmi=spm.inference(0.05);
normality=~spmi.h0reject;

subplot(131);  plot(Y', 'k');  hold on;  title('Data')
subplot(132);  plot(spm.residuals', 'k');  title('Residuals')
subplot(133);  spmi.plot();  title('Normality test')

%--------------------------------------------------------------------------

cond=["Even" "Medium" "High"];
group=["Adultes" "Adolescents" "Enfants" "Jeunes Enfants"];
ind=ones(3,4);
CRPmp=cell(3,4);
for p=1:nbp
    if isempty(PA_CRP{1,p})
        continue
    end
    idg=participants{p,3};
    for c=1:3
        CRPtemp(1,:)=PA_CRP{c,p}{end-idms,4}(paire,:);                      % jambe gauche
        CRPtemp(2,:)=PA_CRP{c+3,p}{end-idms,4}(paire,:);                    % jambe droite
        CRPmp{c,idg}(ind(c,idg),:)=mean(CRPtemp,1);
        ind(c,idg)=ind(c,idg)+1;
    end
end
CRPm=cell(3,4);
for i=1:3
    for j=1:4
        CRPm{i,j}=mean(CRPmp{i,j},1);
        CRPm{i,j}(2,:)=std(CRPmp{i,j},0,1);
    end
end

figure
tl=tiledlayout(3,4);
c=[1 2 3];
cpaire=nchoosek(c,2);
cms=colormap(nebula(3));
cms=brighten(cms,-.2);
for i=1:3
    for g=1:4
        Y0=CRPmp{cpaire(i,1),g};
        Y1=CRPmp{cpaire(i,2),g};
        bc=3;                                                               % Bonferonni correction
rng("shuffle")
        if normality==1
            t=spm1d.stats.ttest_paired(Y0,Y1); 
            ti=t.inference(0.05/bc,'two_tailed',true,'interp',true);
        else
            t=spm1d.stats.nonparam.ttest_paired(Y0,Y1); 
            ti=t.inference(0.05/bc,'two_tailed',true,'iterations',-1,'force_iterations',true);      % iterations : -1 for max iterations
        end
        % For non-parametric tests, a number of iterations smaller than the
        % max may yield qualitatively different results, especially if the
        % effect size is small. Here the max number of iterations is around
        % 30 000, which is fine to run, but if this number was greater, we
        % could fix iterations at 10 000

        nexttile
        hold on
        rect_color(ti,Y0,Y1,0)
        for k=1:2
            plot(CRPm{cpaire(i,k),g}(1,:),'Color',cms(cpaire(i,k),:),'LineWidth',2.5)
            f=fill([1:1:100 100:-1:1],[(CRPm{cpaire(i,k),g}(1,:)+CRPm{cpaire(i,k),g}(2,:)) fliplr((CRPm{cpaire(i,k),g}(1,:)-CRPm{cpaire(i,k),g}(2,:)))],'c');
            f.FaceColor=cms(cpaire(i,k),:);
            f.EdgeColor='none';
            f.FaceAlpha=0.25;
        end
        if idms==1
            ylim([0 180])
            legend(cond(1,cpaire(i,1)),"",cond(1,cpaire(i,2)),"",'Location','southwest')
        else
            ylim([0 60])
            legend(cond(1,cpaire(i,1)),"",cond(1,cpaire(i,2)),"",'Location','northeast')
        end
        box on
        set(gca,'Layer','Top');
        xlabel("% gait cycle")
        if g==1
            if idms==1
                ylabel("MARP (°)")
            else
                ylabel("SD (°)")
            end
        end
        if i==1
            title(group(1,g))
        end
    end
end
if exist('cms_d.mat','file')==2
    load cms_d.mat cms_d
    cms_d=brighten(cms_d,.4);
else
    cms_d=colormap(turbo(300));
    cms_d=brighten(cms_d,.6);
end
colormap(cms_d);
cb=colorbar;
clim([0 3])
cb.Layout.Tile='east';
ylabel(cb,"Cohen's d",'FontSize',13)
if paire==1 && idms==1
    title(tl,"Post-Hoc tests - Knee/Hip MARP",'FontWeight','bold')
elseif paire==1 && idms==0
    title(tl,"Post-Hoc tests - Knee/Hip DP",'FontWeight','bold')
elseif paire==2 && idms==1
    title(tl,"Post-Hoc tests - Ankle/Knee MARP",'FontWeight','bold')
elseif paire==2 && idms==0
    title(tl,"Post-Hoc tests - Ankle/Knee DP",'FontWeight','bold')
end
tl.Padding = 'compact'; tl.TileSpacing = 'compact';

figure
tl=tiledlayout(3,6);
g=[1 2 3 4];
gpaire=nchoosek(g,2);
cms=colormap(turbo(4));
for c=1:3
    for i=1:6
        Y0=CRPmp{c,gpaire(i,1)};
        Y1=CRPmp{c,gpaire(i,1)};
        bc=6;                                                               % Bonferonni correction

        if normality==1
            t=spm1d.stats.ttest_paired(Y0,Y1); 
            ti=t.inference(0.05/bc,'two_tailed',true,'interp',true);
        else
            t=spm1d.stats.nonparam.ttest_paired(Y0,Y1); 
            ti=t.inference(0.05/bc,'two_tailed',true,'iterations',-1,'force_iterations',true);
        end

        nexttile
        hold on
        rect_color(ti,Y0,Y1,0)
        for k=1:2
            plot(CRPm{c,gpaire(i,k)}(1,:),'Color',cms(gpaire(i,k),:),'LineWidth',2.5)
            f=fill([1:1:100 100:-1:1],[(CRPm{c,gpaire(i,k)}(1,:)+CRPm{c,gpaire(i,k)}(2,:)) fliplr((CRPm{c,gpaire(i,k)}(1,:)-CRPm{c,gpaire(i,k)}(2,:)))],'c');
            f.FaceColor=cms(gpaire(i,k),:);
            f.EdgeColor='none';
            f.FaceAlpha=0.2;
        end
        if idms==1
            ylim([0 180])
            legend(group(1,gpaire(i,1)),"",group(1,gpaire(i,2)),"",'Location','southwest')
        else
            ylim([0 60])
            legend(group(1,gpaire(i,1)),"",group(1,gpaire(i,2)),"",'Location','northeast')
        end
        box on
        set(gca,'Layer','Top');
        xlabel("% gait cycle")
        if i==1
            ylabel(cond(1,c),'FontWeight','bold')
        end
        ysecondarylabel("(°)")
    end
end
if paire==1 && idms==1
    title(tl,"Post-Hoc tests - Knee/Hip MARP",'FontWeight','bold')
elseif paire==1 && idms==0
    title(tl,"Post-Hoc tests - Knee/Hip DP",'FontWeight','bold')
elseif paire==2 && idms==1
    title(tl,"Post-Hoc tests - Ankle/Knee MARP",'FontWeight','bold')
elseif paire==2 && idms==0
    title(tl,"Post-Hoc tests - Ankle/Knee DP",'FontWeight','bold')
end
tl.Padding = 'compact'; tl.TileSpacing = 'compact';

%% Plan de covariation

p=65;
c=3;

load cov_pla.mat
at=cov_pla{c,p}{3,1}(1,:);
as=cov_pla{c,p}{3,1}(2,:);
af=cov_pla{c,p}{3,1}(3,:);

figure
cms=colormap(nebula(100));
scatter3(at,as,af,50,cms,'filled')
xlabel("\theta thigh (°)"); ylabel("\theta shank (°)"); zlabel("\theta foot (°)")
hold on

v1=cov_pla{c,p}{4,1}(:,1);
v2=cov_pla{c,p}{4,1}(:,2);
p0=[mean(at) mean(as) mean(af)];

s_range=-ceil(abs(min(cov_pla{c,p}{5,1}(:,1)))/5)*5-5:5:ceil(max(cov_pla{c,p}{5,1}(:,1))/5)*5+5;
t_range=-ceil(abs(min(cov_pla{c,p}{5,1}(:,2)))/5)*5-5:5:ceil(max(cov_pla{c,p}{5,1}(:,2))/5)*5+5;
[s,t]=meshgrid(s_range,t_range);
X=p0(1)+s*v1(1)+t*v2(1);
Y=p0(2)+s*v1(2)+t*v2(2);
Z=p0(3)+s*v1(3)+t*v2(3);
surf(X,Y,Z,"FaceColor","none","EdgeColor",[.3 .3 .3]);
xlim([min(min(X,[],"all"),min(at)) max(max(X,[],"all"),max(at))])
ylim([min(min(Y,[],"all"),min(as)) max(max(Y,[],"all"),max(as))])
zlim([min(min(Z,[],"all"),min(af)) max(max(Z,[],"all"),max(af))])
% ev1=mean(abs(cov_pla{c,p}{5,1}(:,1)));
% ev2=mean(abs(cov_pla{c,p}{5,1}(:,2)));
q=quiver3(p0(1),p0(2),p0(3),v1(1)*15,v1(2)*15,v1(3)*15,'r','linewidth',2);
q.MaxHeadSize=.5;
q=quiver3(p0(1),p0(2),p0(3),v2(1)*15,v2(2)*15,v2(3)*15,'b','linewidth',2);
q.MaxHeadSize=.5;
title("Covariation plane of elevation angles")
legend("","","PC1","PC2",Location="northeast")
box on

%% Plan de covariation - u3t - 3 conditions

p=4;

load cov_pla.mat
figure
tl=tiledlayout(1,3);
cms=colormap(nebula(100));
cond=["Even" "Medium" "High"];

for c=1:3
    nexttile
    at=cov_pla{c,p}{3,1}(1,:);
    as=cov_pla{c,p}{3,1}(2,:);
    af=cov_pla{c,p}{3,1}(3,:);
    scatter3(at,as,af,50,cms,'filled')
    xlabel("\theta thigh (°)"); ylabel("\theta shank (°)"); zlabel("\theta foot (°)")
    hold on
    v1=cov_pla{c,p}{4,1}(:,1);
    v2=cov_pla{c,p}{4,1}(:,2);
    p0=[mean(at) mean(as) mean(af)];
    s_range=-ceil(abs(min(cov_pla{c,p}{5,1}(:,1)))/5)*5-5:5:ceil(max(cov_pla{c,p}{5,1}(:,1))/5)*5+5;
    t_range=-ceil(abs(min(cov_pla{c,p}{5,1}(:,2)))/5)*5-5:5:ceil(max(cov_pla{c,p}{5,1}(:,2))/5)*5+5;
    [s,t]=meshgrid(s_range,t_range);
    X=p0(1)+s*v1(1)+t*v2(1);
    Y=p0(2)+s*v1(2)+t*v2(2);
    Z=p0(3)+s*v1(3)+t*v2(3);
    surf(X,Y,Z,"FaceColor","none","EdgeColor",[.3 .3 .3]);
    xlim([min(min(X,[],"all"),min(at)) max(max(X,[],"all"),max(at))])
    ylim([min(min(Y,[],"all"),min(as)) max(max(Y,[],"all"),max(as))])
    zlim([min(min(Z,[],"all"),min(af)) max(max(Z,[],"all"),max(af))])
    title(cond(c))
    view(-90,0)
    box on
end
title(tl,"Orientation of the covariation plane of elevation angles")
tl.Padding = 'compact'; tl.TileSpacing = 'compact';