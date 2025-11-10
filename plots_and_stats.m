
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
        CRPtemp(1,:)=PA_CRP{c,p}{end-idms,4}(paire,:);                         % jambe gauche
        CRPtemp(2,:)=PA_CRP{c+3,p}{end-idms,4}(paire,:);                       % jambe droite
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
cms=colormap(parula(3));
for i=1:3
    for g=1:4
        Y0=CRPmp{cpaire(i,1),g};
        Y1=CRPmp{cpaire(i,2),g};
        bc=3;                                                               % bonferonni correction
        t=spm1d.stats.ttest_paired(Y0,Y1); 
        ti=t.inference(0.05/bc,'two_tailed',true,'interp',true);            % CHECKER PARAMETRES

        nexttile
        hold on
        nbclust=size(ti.clusters,2);
        for j=1:nbclust
            rectangle('Position',[ti.clusters{1,j}.endpoints(1,1) 0 ti.clusters{1,j}.endpoints(1,2)-ti.clusters{1,j}.endpoints(1,1) 180],'EdgeColor','none','FaceColor',[.9 .9 .9])
        end
        for k=1:2
            plot(CRPm{cpaire(i,k),g}(1,:),'Color',cms(cpaire(i,k),:),'LineWidth',2.5)
            f=fill([1:1:100 100:-1:1],[(CRPm{cpaire(i,k),g}(1,:)+CRPm{cpaire(i,k),g}(2,:)) fliplr((CRPm{cpaire(i,k),g}(1,:)-CRPm{cpaire(i,k),g}(2,:)))],'c');
            f.FaceColor=cms(cpaire(i,k),:);
            f.EdgeColor='none';
            f.FaceAlpha=0.2;
        end
        ylim([0 180])
        legend(cond(1,cpaire(i,1)),"",cond(1,cpaire(i,2)),"",'Location','southwest')
        xlabel("% gait cycle")
        if g==1
            ylabel("MARP (°)")
        end
        if i==1
            title(group(1,g))
        end
    end
end
if paire==1
    title(tl,"Post-Hoc tests - Knee/Hip MARP",'FontWeight','bold')
else
    title(tl,"Post-Hoc tests - Ankle/Knee MARP",'FontWeight','bold')
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
        bc=6;                                                               % bonferonni correction
        t=spm1d.stats.ttest2(Y0,Y1); 
        ti=t.inference(0.05/bc,'two_tailed',true,'interp',true);            % CHECKER PARAMETRES

        nexttile
        hold on
        nbclust=size(ti.clusters,2);
        for j=1:nbclust
            rectangle('Position',[ti.clusters{1,j}.endpoints(1,1) 0 ti.clusters{1,j}.endpoints(1,2)-ti.clusters{1,j}.endpoints(1,1) 180],'EdgeColor','none','FaceColor',[.9 .9 .9])
        end
        for k=1:2
            plot(CRPm{c,gpaire(i,k)}(1,:),'Color',cms(gpaire(i,k),:),'LineWidth',2.5)
            f=fill([1:1:100 100:-1:1],[(CRPm{c,gpaire(i,k)}(1,:)+CRPm{c,gpaire(i,k)}(2,:)) fliplr((CRPm{c,gpaire(i,k)}(1,:)-CRPm{c,gpaire(i,k)}(2,:)))],'c');
            f.FaceColor=cms(gpaire(i,k),:);
            f.EdgeColor='none';
            f.FaceAlpha=0.2;
        end
        ylim([0 180])
        legend(group(1,gpaire(i,1)),"",group(1,gpaire(i,2)),"",'Location','southwest')
        xlabel("% gait cycle")
        if i==1
            ylabel(cond(1,c),'FontWeight','bold')
        end
        ysecondarylabel("(°)")
    end
end
if paire==1
    title(tl,"Post-Hoc tests - Knee/Hip MARP",'FontWeight','bold')
else
    title(tl,"Post-Hoc tests - Ankle/Knee MARP",'FontWeight','bold')
end
tl.Padding = 'compact'; tl.TileSpacing = 'compact';

%% Plan de covariation

p=4;
c=1;

load cov_pla.mat
at=cov_pla{c,p}{3,1}(1,:);
as=cov_pla{c,p}{3,1}(2,:);
af=cov_pla{c,p}{3,1}(3,:);

figure
cms=colormap(nebula(100));
scatter3(at,as,af,50,cms,'filled')
xlabel("\theta thigh"); ylabel("\theta shank"); zlabel("\theta foot")
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
ev1=mean(abs(cov_pla{c,p}{5,1}(:,1)));
ev2=mean(abs(cov_pla{c,p}{5,1}(:,2)));
q=quiver3(p0(1),p0(2),p0(3),v1(1)*ev1,v1(2)*ev1,v1(3)*ev1,'r','linewidth',2);
q.MaxHeadSize=.5;
q=quiver3(p0(1),p0(2),p0(3),v2(1)*ev2,v2(2)*ev2,v2(3)*ev2,'b','linewidth',2);
q.MaxHeadSize=.5;
