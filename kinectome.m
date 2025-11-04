% Ce script utilise le cell array Acc fourni par mark_extraction pour créer
% un kinectome de chaque essai pour chaque participant. Ces kinectomes
% ainsi que les moyennes et écarts-types de chaque paire part/cond sont
% stockés dans le cell array Kinect. Le 'weighted degree" de chaque node 
% est calculé sur les kinectomes moyens et stockés dans le cell array WD.

%%
clc
clear
%--------------------------------------------------------------------------
load Acc.mat
nbp=size(Acc,2);
Kinect=cell(3,nbp);
WD=cell(3,nbp);

for p=1:nbp
    if isempty(Acc{1,p})
        continue
    end
    for c=1:3
        nbc=size(Acc{c,p},2)/100;
        for cy=1:nbc
            atemp=Acc{c,p}(:,1+(cy-1)*100:100*cy);
            ahead=atemp(1:8,:);                                             % Acc. medlat et antpos des 4 marqueurs de tête
            ahead(all(ahead==0,2),:)=[];                                    % Marqueur manquant pour ce cycle
            ahml=mean(ahead(1:2:end,:),1);
            ahap=mean(ahead(2:2:end,:),1);                                  % Acc. medlat et antpos moyennes de la tête
            atemp(1:8,:)=[];
            a_glob=[ahml;ahap;atemp];
            a_medlat=a_glob(1:2:end,:);
            a_antpos=a_glob(2:2:end,:);
            k_medlat=corrcoef(a_medlat');                                   % Kinectome par corrélation de Pearson
            k_antpos=corrcoef(a_antpos');
            Kinect{c,p}{cy,1}=k_medlat;
            Kinect{c,p}{cy,2}=k_antpos;

            if ~exist('Kml','var') && ~any(isnan(k_medlat(:)))              % Le Kinectome moyen ne prend pas en compte les cycles avec un marqueur manquant
                Kml=k_medlat;
                Kap=k_antpos;  
            elseif ~any(isnan(k_medlat(:)))
                Kap=cat(3,Kap,k_antpos);
                Kml=cat(3,Kml,k_medlat);
            end
        end
        Kinect{c,p}{nbc+1,2}=mean(Kml,3);                                   % Moyenne des kinectomes en avant-dernière ligne
        Kinect{c,p}{nbc+2,2}=std(Kml,0,3);                                  % Ecarts types en dernière ligne
        Kinect{c,p}{nbc+1,1}=mean(Kap,3);
        Kinect{c,p}{nbc+2,1}=std(Kap,0,3);
        clear Kap Kml
        nbm=size(k_antpos,1);
        weig_deg=[nbm,2];                                                   % weighted degree of nodes
        for node=1:nbm
            weig_deg(node,1)=(sum(abs(Kinect{c,p}{end-1,1}(node,:)))-1)/nbm;
            weig_deg(node,2)=(sum(abs(Kinect{c,p}{end-1,2}(node,:)))-1)/nbm;
        end
        WD{c,p}=weig_deg;
    end
end

save Kinect.mat Kinect
save WD.mat WD