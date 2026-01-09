% Ce script identifie les cycles de marche entiers de chaque essai puis
% extrait les angles en flex/ext hanche, genou et cheville. Pour chaque
% participant et chaque condition, les angles  et les %TO sont regroupés
% dans le cell array 'K_crp'
%
%%
clc;
clear;
close all;
% -------------------------------------------------------------------------

addpath('.\btk');
nbp=70;                                                                     % Nombre de participants
cond={'Plat' 'Medium' 'High'};
nbe=10;                                                                     % Nombre d'essais
ang={'LHipAngles' 'LKneeAngles' 'LAnkleAngles';...
     'RHipAngles' 'RKneeAngles' 'RAnkleAngles'};

% -------------------------------------------------------------------------
if exist('K_crp.mat','file')==2
    load K_crp.mat
    nbpk=size(K_crp,2);
    if nbpk<nbp
        K_crp=[K_crp cell(6,nbp-nbpk)];
    end
else
    K_crp=cell(6,nbp);
    % 6 Lignes : Gauche ('Plat' 'Medium' 'High'), Droite ('Plat' 'Medium' 'High')
end

[bf,af]=butter(4,6/(100/2),'low');

for p=1:nbp
    part=sprintf('CTL_%02d',p);
    disp(['Processing participant: ' part]);
    temp=[part '_Plat_01.c3d'];
    if ~isempty(K_crp{1,p})
        continue
    elseif ~exist(temp,'file')
        continue
    end
    for c=1:length(cond)
        disp(['Condition: ' cond{c}]);
        for j=1:2                                                           % Jambe g/d
            ind=1;
            for e=1:nbe
                ess=sprintf('%02d',e);
                file=[part '_' cond{c} '_' ess '.c3d'];
                if ~exist(file,'file')
                    continue
                end
                data=btkReadAcquisition(file);
                angles=btkGetAngles(data);
                for a=1:length(ang)
                    angles.(ang{j,a})=filtfilt(bf,af,angles.(ang{j,a}));
                end
                events=btkGetEvents(data);
                start=btkGetFirstFrame(data);
                if j==1
                    HS=round(events.Left_Foot_Strike*100-start);            % Heel strikes
                    TO=round(events.Left_Foot_Off*100-start);               % Toe offs
                else
                    HS=round(events.Right_Foot_Strike*100-start);
                    TO=round(events.Right_Foot_Off*100-start);
                end
                HS(HS<=0)=1;                                                % Au cas où HS1 coïncide avec la première frame ou a une frame d'avance
                nbc=length(HS)-1;                                           % Nombre de cycles entiers
                for cy=1:nbc
                    ma=[];
                    for a=1:length(ang)
                        va=angles.(ang{j,a})(HS(cy):HS(cy+1),1);
                        ma=[ma;va'];
                    end
                    TOcy=TO(TO>HS(cy) & TO<HS(cy+1));
                    if isempty(TOcy)
                        TOperc=0;
                    else
                        TOperc=((TOcy-HS(cy))/(HS(cy+1)-HS(cy)))*100;       % Pourcentage Toe Off
                    end
                    % disp(['TO : ' num2str(TOperc)])
                    K_crp{c+(j-1)*3,p}{ind,1}=ma;
                    K_crp{c+(j-1)*3,p}{ind,2}=TOperc;
                    ind=ind+1;
                end
            end
        end
    end
end

save K_crp.mat K_crp