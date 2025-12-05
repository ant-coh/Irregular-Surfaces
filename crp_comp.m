% Ce script calcule les angles de phases des trois articulations puis la
% CRP hanche-genou et genou-cheville pour chaque cycle de marche. La MARP
% et l'écart-type des crp sont ensuite calculés pour chaque p et c.

%%
clc;
clear;
% -------------------------------------------------------------------------
addpath('Functions')
load K_crp.mat
nbp=length(K_crp);
PA_CRP=cell(6,nbp);
% Lignes : Gauche ('Plat' 'Medium' 'High'), Droite ('Plat' 'Medium' 'High')
% Cellules : Colonnes PA PAn CRP CRPm (MARP et std en dernières lignes)

for p=1:nbp
    if isempty(K_crp{1,p})
        continue
    end
    for c=1:3
        for j=1:2
            nbc=length(K_crp{c+(j-1)*3,p});
            for cy=1:nbc
                for a=1:3
                    ang=K_crp{c+(j-1)*3,p}{cy,1}(a,:);
                    pad=round(0.1*length(ang));
                    x=1:1:length(ang);
                    xq=1-pad:1:length(ang)+pad;
                    pad_ang=spline(x,ang,xq);
                    PA=Phase_Angle(pad_ang');
                    PA_CRP{c+(j-1)*3,p}{cy,1}(a,:)=(PA(pad+1:end-pad,1))';
                end
                PAtemp=PA_CRP{c+(j-1)*3,p}{cy,1};
                PAn=interp1(linspace(1,size(PAtemp,2),size(PAtemp,2)),PAtemp',linspace(1,size(PAtemp,2),101));
                PA_CRP{c+(j-1)*3,p}{cy,2}=PAn';
                if ~exist('PAm','var')
                    PAm=PAn';
                else
                    PAm=cat(3,PAm,PAn');
                end    

                HKcrp=CRP(PA_CRP{c+(j-1)*3,p}{cy,1}(2,:),PA_CRP{c+(j-1)*3,p}{cy,1}(1,:));
                KAcrp=CRP(PA_CRP{c+(j-1)*3,p}{cy,1}(3,:),PA_CRP{c+(j-1)*3,p}{cy,1}(2,:));
                PA_CRP{c+(j-1)*3,p}{cy,3}(1,:)=HKcrp;
                PA_CRP{c+(j-1)*3,p}{cy,3}(2,:)=KAcrp;
                CRPtemp=PA_CRP{c+(j-1)*3,p}{cy,3};
                CRPn=interp1(linspace(1,size(CRPtemp,2),size(CRPtemp,2)),CRPtemp',linspace(1,size(CRPtemp,2),101));
                PA_CRP{c+(j-1)*3,p}{cy,4}=CRPn';
                if ~exist('CRPm','var')
                    CRPm=CRPn';
                else
                    CRPm=cat(3,CRPm,CRPn');
                end

            end
            PA_CRP{c+(j-1)*3,p}{nbc+1,2}=mean(PAm,3);
            PA_CRP{c+(j-1)*3,p}{nbc+2,2}=std(PAm,0,3);
            PA_CRP{c+(j-1)*3,p}{nbc+1,4}=mean(CRPm,3);
            PA_CRP{c+(j-1)*3,p}{nbc+2,4}=std(CRPm,0,3);
            clear PAm CRPm
        end
    end
end

save PA_CRP.mat PA_CRP