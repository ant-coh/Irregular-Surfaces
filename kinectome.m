% Ce script utilise le cell array Acc fourni par mark_extraction pour créer
% un kinectome de chaque essai pour chaque participant. Ces kinectomes
% ainsi que les moyennes et écarts-types de chaque paire part/cond sont
% stockés dans le cell array Kinect. La corr. abs. moyenne de chaque node 
% est calculée sur les kinectomes moyens et stockés dans le cell array MAC.

%%
clc
clear
%--------------------------------------------------------------------------
load Acc.mat
nbp=size(Acc,2);
Kinect=cell(3,nbp);
MAC=cell(3,nbp);
sil=cell(3,nbp);

S=[1 1 1 1 2 2 1 3 3 1 1 3 3 3 3 3 2 2 2 2 2];
[S_sort,I]=sort(S);
nbclust=max(S_sort);

for p=1:nbp
    if isempty(Acc{1,p})
        continue
    end
    for c=1:3
        nbc=size(Acc{c,p},2)/101;
        for cy=1:nbc
            atemp=Acc{c,p}(:,1+(cy-1)*101:101*cy);
            ahead=atemp(1:8,:);                                             % Acc. medlat et antpos des 4 marqueurs de tête
            ahead(all(ahead==0,2),:)=[];                                    % Marqueur manquant pour ce cycle
            ahml=mean(ahead(1:2:end,:),1);
            ahap=mean(ahead(2:2:end,:),1);                                  % Acc. medlat et antpos moyennes de la tête
            atemp(1:8,:)=[];
            a_glob=[ahml;ahap;atemp];

            a_glob(53:54,:)=[];                                             % Suppression LHEE et RHEE
            a_glob(41:42,:)=[];
            a_glob(29:32,:)=[];                                             % Suppression LPSI RPSI
            a_glob(7:12,:)=[];                                              % Suppression CLAV, STRN, RBAK

            a_medlat=a_glob(1:2:end,:);
            a_antpos=a_glob(2:2:end,:);
            
            k_medlat=corrcoef(a_medlat');                                   % Kinectome par corrélation de Pearson
            k_antpos=corrcoef(a_antpos');
            Kinect{c,p}{cy,2}=k_medlat;
            Kinect{c,p}{cy,1}=k_antpos;

            if ~exist('Kml','var') && ~any(isnan(k_medlat(:)))              % Le Kinectome moyen ne prend pas en compte les cycles avec un marqueur manquant
                Kml=k_medlat;
                Kap=k_antpos;
            elseif ~any(isnan(k_medlat(:)))
                Kap=cat(3,Kap,k_antpos);
                Kml=cat(3,Kml,k_medlat);
            end
        end
        Kinect{c,p}{nbc+1,2}=mean(Kml,3);                                   % Moyenne des kinectomes en avant-dernière ligne
        Kinect{c,p}{nbc+2,2}=std(Kml,0,3);                                  % Ecarts types en dernière ligne (peu représentatifs car les cycles se chevauchent)
        Kinect{c,p}{nbc+1,1}=mean(Kap,3);
        Kinect{c,p}{nbc+2,1}=std(Kap,0,3);
        clear Kap Kml
        nbm=size(k_antpos,1);
        mac=[nbm,2];                                                        % Mean absolute correlation
        silhou=zeros(1,nbm);

        for i=1:nbm
            mac(i,1)=(sum(abs(Kinect{c,p}{end-1,1}(i,:)))-1)/nbm;
            mac(i,2)=(sum(abs(Kinect{c,p}{end-1,2}(i,:)))-1)/nbm;
            
            Kin=Kinect{c,p}{end-1,1}(I,I);
            clust_i=S_sort(i);
            ind=1;
            temp=S_sort;
            dist=Kin(i,:);
            dist(i)=[];
            temp(i)=[];
            b=zeros(1,2);
            for cl=1:nbclust
                if cl==clust_i
                    v=dist(temp==clust_i);
                    a=1-mean(v);
                else
                    v=dist(temp==cl);
                    b(ind)=1-mean(v);
                    ind=ind+1;
                end
            end
            silhou(i)=(min(b)-a)/max(a,min(b));
        end
        MAC{c,p}=mac;
        sil{c,p}=mean(silhou,2);
    end
end

save Kinect.mat Kinect
save MAC.mat MAC
save sil.mat sil