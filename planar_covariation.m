% This script computes elevations angles for each cycle (left and right) of
% each condition for each participant. PCA is then performed on the means
% of these angles for each participant and condition.

%%
clear
clc
%--------------------------------------------------------------------------

load K.mat
load K_crp.mat
nbp=size(K,2);
cov_pla=cell(3,nbp);
addpath('Functions\')

for p=1:nbp
    if isempty(K{1,p})
        continue
    end
    for c=1:3
        for j=1:2
            nbc=length(K{c+3*(j-1),p})/100;
            for cy=1:nbc
                Ktemp=K{c+3*(j-1),p}(:,(1+(cy-1)*100):100*cy);
                ap=Ktemp(16,:);                                             % Pelvis forward tilt / global axes
                ah=Ktemp(1,:);                                              % Hip f/e
                ak=Ktemp(4,:);                                              % Knee f/e
                aa=Ktemp(7,:);                                              % Ankle f/e
                at=ah-ap;                                                   % Thigh elevation angle
                as=at-ak;                                                   % Shank
                af=90+as+aa;                                                % Foot
                if cy==1
                    mat=at;
                    mas=as;
                    maf=af;
                else
                    mat=[mat;at];
                    mas=[mas;as];
                    maf=[maf;af];
                end
            end
            atm=mean(mat,1);                                                % Mean elevation angles
            asm=mean(mas,1);
            afm=mean(maf,1);
            cov_pla{c,p}{j,1}=[reshape(mat.',1,[]);reshape(mas.',1,[]);reshape(maf.',1,[])];
            cov_pla{c,p}{j+2,1}=[atm;asm;afm];
            for cy=1:nbc
                PAtemp=zeros(3,100);
                for a=1:3
                    ang=cov_pla{c,p}{j,1}(a,1+(cy-1)*100:100*cy);
                    pad=round(0.1*length(ang));
                    x=1:1:length(ang);
                    xq=1-pad:1:length(ang)+pad;
                    pad_ang=spline(x,ang,xq);
                    PA=Phase_Angle(pad_ang');
                    PAtemp(a,:)=(PA(pad+1:end-pad,1))';
                end
                crp_ts=CRP(PAtemp(2,:),PAtemp(1,:));
                crp_sf=CRP(PAtemp(3,:),PAtemp(2,:));
                if cy==1
                    mpt=PAtemp(1,:);
                    mps=PAtemp(2,:);
                    mpf=PAtemp(3,:);
                    m_crpts=crp_ts;
                    m_crpsf=crp_sf;
                else
                    mpt=[mpt;PAtemp(1,:)];
                    mps=[mps;PAtemp(2,:)];
                    mpf=[mpf;PAtemp(3,:)];
                    m_crpts=[m_crpts;crp_ts];
                    m_crpsf=[m_crpsf;crp_sf];
                end
            end
            ptm=mean(mpt,1);                                                % Mean phase angles
            psm=mean(mps,1);
            pfm=mean(mpf,1);
            cov_pla{c,p}{j,2}=[reshape(mpt.',1,[]);reshape(mps.',1,[]);reshape(mpf.',1,[])];
            cov_pla{c,p}{j+2,2}=[ptm;psm;pfm];
        end
        temp=cat(3,cov_pla{c,p}{3,1},cov_pla{c,p}{4,1});
        mtemp=mean(temp,3);
        cov_pla{c,p}{3,1}=[];
        cov_pla{c,p}{4,1}=[];
        cov_pla{c,p}{3,1}=mtemp;
        [PC,score,latent,~,Vpc]=pca([mtemp(1,:)' mtemp(2,:)' mtemp(3,:)']); % PCA
        cov_pla{c,p}{4,1}=PC;                                               % Eigenvectors
        cov_pla{c,p}{5,1}=score;
        cov_pla{c,p}{6,1}=latent;                                           % Eigenvalues
        cov_pla{c,p}{7,1}=Vpc;                                              % Variance explained
        cov_pla{c,p}{8,1}=Vpc(1)+Vpc(2);                                    % Planarity index
        temp=cat(3,cov_pla{c,p}{3,2},cov_pla{c,p}{4,2});
        mtemp=mean(temp,3);
        cov_pla{c,p}{3,2}=[];
        cov_pla{c,p}{4,2}=[];
        cov_pla{c,p}{3,2}=mtemp;
        cov_pla{c,p}{4,2}(1,:)=mean(m_crpts,1);
        cov_pla{c,p}{4,2}(2,:)=mean(m_crpsf,1);
        cov_pla{c,p}{5,2}=(mean(cell2mat(K_crp{c,p}(:,2)),1)+mean(cell2mat(K_crp{c+3,p}(:,2)),1))/2; % Mean toe off percentage
        for i=1:2
            cov_pla{c,p}{6,2}(i,1)=mean(cov_pla{c,p}{4,2}(i,:),2);
            TO=(round(cov_pla{c,p}{5,2}));                                  % A vérifier
            cov_pla{c,p}{7,2}(i,1)=mean(cov_pla{c,p}{4,2}(i,1:TO),2);
            cov_pla{c,p}{8,2}(i,1)=mean(cov_pla{c,p}{4,2}(i,TO:100),2);
        end
    end
end

save cov_pla.mat cov_pla