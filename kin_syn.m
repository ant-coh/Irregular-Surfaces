
clear
clc
%--------------------------------------------------------------------------

load K.mat
nbp=size(K,2);
W=cell(6,nbp);
C=cell(6,nbp);
V=cell(6,nbp);
Km=cell(6,nbp);
Ks=cell(6,nbp);

for p=1:nbp
    if isempty(K{1,p})
        continue
    end
    k=K(1:6,p);
    [nba,~]=size(k{1,1});
    for j=1:2                                                               % Normalisation par la valeur absolue de l'angle max
        for i=1:nba                                                         % atteint au cours des trois essais
            amax=max([max(abs(k{1+3*(j-1)}(i,:))) max(abs(k{2+3*(j-1)}(i,:))) max(abs(k{3+3*(j-1)}(i,:)))]);
            for c=1:3
                k{c+3*(j-1)}(i,:)=k{c+3*(j-1)}(i,:)/amax;
            end
        end
        for c=1:3
            for cy=1:length(k{c+3*(j-1),1})/101
                Ktemp=k{c+3*(j-1),1}(:,(1+(cy-1)*101):101*cy);
                if ~exist('Ktm','var')
                    Ktm=Ktemp;
                else
                    Ktm=cat(3,Ktm,Ktemp);
                end
            end
            Km{c+3*(j-1),p}=mean(Ktm,3);                                    % Mean cycles
            Ks{c+3*(j-1),p}=std(Ktm,0,3);                                   % Standard deviation
            clear Ktm
        end
    end

%--------------------------------------------------------------------------
    ang_sup=[10 11 12 13 16 18];                                            % Suppression angles non-pertinents
%    ang_sup=[10 11 12 13 14 15 16 17 18];
    for j=1:2
        for c=1:3
            for i=1:length(ang_sup)
                k{c+3*(j-1)}(ang_sup(i)-i+1,:) = [];
            end
            nbc=size(k{c+3*(j-1)},2)/101;
            for cy=1:nbc
                [Wtemp,Ctemp,~,~,Vtemp]=pca((k{c+3*(j-1)}(:,1+(cy-1)*101:101*cy))');
                for i=1:nba-length(ang_sup)
                    if cy==1 && c~=1 && min(corrcoef(Ctemp(:,i),C{1+3*(j-1),p}{1,end-1}(i,:)),[],"all")<0
                        Ctemp(:,i)=-Ctemp(:,i);
                        Wtemp(:,i)=-Wtemp(:,i);
                    elseif cy~=1 && min(corrcoef(Ctemp(:,i),Ccomp(:,i)),[],"all")<0
                        Ctemp(:,i)=-Ctemp(:,i);
                        Wtemp(:,i)=-Wtemp(:,i);
                    end
                end
                Ccomp=Ctemp;
                W{c+3*(j-1),p}{1,cy}=Wtemp;
                C{c+3*(j-1),p}{1,cy}=Ctemp';
                V{c+3*(j-1),p}{1,cy}=Vtemp;
                if cy==1
                    Wm=Wtemp;
                    Cm=Ctemp';
                    Vm=Vtemp;
                else
                    Wm=cat(3,Wm,Wtemp);
                    Cm=cat(3,Cm,Ctemp');
                    Vm=cat(3,Vm,Vtemp);
                end
            end
            W{c+3*(j-1),p}{1,nbc+1}=mean(Wm,3);
            W{c+3*(j-1),p}{1,nbc+2}=std(Wm,0,3);
            C{c+3*(j-1),p}{1,nbc+1}=mean(Cm,3);
            C{c+3*(j-1),p}{1,nbc+2}=std(Cm,0,3);
            V{c+3*(j-1),p}{1,nbc+1}=mean(Vm,3);
            V{c+3*(j-1),p}{1,nbc+2}=std(Vm,0,3);
        end
    end
end

%%