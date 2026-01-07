% Ce script identifie les cycles de marche entiers de chaque essai puis
% extrait les angles articulaires pertinents. Pour chaque participant et
% chaque condition, les angles sont normalisés sur un cycle (100%) puis
% concaténés. Toutes les données sont regroupées dans le cell array K.
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
ang={'LHipAngles' 'LKneeAngles' 'LAnkleAngles' 'LThoraxAngles' 'LSpineAngles' 'LPelvisAngles';...
     'RHipAngles' 'RKneeAngles' 'RAnkleAngles' 'RThoraxAngles' 'RSpineAngles' 'RPelvisAngles'};
% 3 angles pour chaque articulation : Flex/Ext, Add/Abs, Rot
%                                 ou Fwd tilt, Lat tilt, Rot

% -------------------------------------------------------------------------
if exist('K.mat','file')==2
    load K.mat
    nbpk=size(K,2);
    if nbpk<nbp
        K=[K cell(6,nbp-nbpk)];
    end
else
    K=cell(6,nbp);
    % 6 Lignes : Gauche ('Plat' 'Medium' 'High'), Droite ('Plat' 'Medium' 'High')
end

[b,a]=butter(4,6/(100/2),'low');                                            % Filtre passe-bas angles articulaires

for p=2:nbp
    part=sprintf('CTL_%02d',p);
    disp(['Processing participant: ' part]);
    temp=[part '_Plat_01.c3d'];
    if ~isempty(K{1,p})
        continue
    elseif ~exist(temp,'file')
        continue
    end
    for c=1:length(cond)
        for j=1:2                                                           % Jambe g/d
            mk=[];
            for e=1:nbe
                ess=sprintf('%02d',e);
                file=[part '_' cond{c} '_' ess '.c3d'];
                if ~exist(file,'file')
                    continue
                end
                data=btkReadAcquisition(file);
                angles=btkGetAngles(data);
                for a=1:length(ang)
                    angles.(ang{j,a})=filtfilt(b,a,angles.(ang{j,a}));
                end
                events=btkGetEvents(data);
                start=btkGetFirstFrame(data);
                if j==1
                    HS=round(events.Left_Foot_Strike*100-start);            % Heel strikes
                else
                    HS=round(events.Right_Foot_Strike*100-start);
                end
                HS(HS<=0)=1;                                                % Au cas où HS1 coïncide avec la première frame ou a une frame d'avance
                nbc=length(HS)-1;                                           % Nombre de cycles entiers
                for cy=1:nbc
                    ma=[];
                    for a=1:length(ang)
                        va=angles.(ang{j,a})(HS(cy):HS(cy+1),:);
                        ma=[ma va];
                    end                                                     % Normalisation des cycles
                    mn=interp1(linspace(1,size(ma,1),size(ma,1)),ma,linspace(1,size(ma,1),101));
                    mk=[mk mn'];                                            % Concaténation
                end
            end
            K{c+(j-1)*3,p}=mk;
        end
        disp(['Condition: ' cond{c}]);
    end
end
save K.mat K