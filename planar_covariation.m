% This script computes elevations angles for each cycle (left and right) of
% each condition for each participant. PCA is then performed on the means
% of these angles for each participant and condition.

%%
clear
clc
%--------------------------------------------------------------------------

load K.mat
nbp=size(K,2);
cov_pla=cell(3,nbp);

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
                ah=Ktemp(1,:);                                              % Hip
                ak=Ktemp(4,:);                                              % Knee
                aa=Ktemp(7,:);                                              % Ankle
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
            cov_pla{c,p}{j,1}=[reshape(mat,1,[]);reshape(mas,1,[]);reshape(maf,1,[])];
            cov_pla{c,p}{j+2,1}=[atm;asm;afm];
        end
        temp=cat(3,cov_pla{c,p}{3,1},cov_pla{c,p}{4,1});
        mtemp=mean(temp,3);
        cov_pla{c,p}{3,1}=[];
        cov_pla{c,p}{4,1}=[];
        cov_pla{c,p}{3,1}=mtemp;
        [PC,score,~,~,Vpc]=pca([mtemp(1,:)' mtemp(2,:)' mtemp(3,:)']);      % PCA
        cov_pla{c,p}{4,1}=PC;                                               % Eigenvectors
        cov_pla{c,p}{5,1}=score;                                            % Eigenvalues
        cov_pla{c,p}{6,1}=Vpc;                                              % Variance explained
        cov_pla{c,p}{7,1}=Vpc(1)+Vpc(2);                                    % Planarity index
    end
end

save cov_pla.mat cov_pla