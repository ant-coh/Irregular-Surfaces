function rect_color(SPMi,Y0,Y1,varargin)

% Plots rectangles between endpoints of clusters identified by spm1d.
% These rectangles are colored based on the mean Cohen's d of the cluster.

if nargin==3
    thres=0;                                                                % Minimum length of displayed clusters
    maxd=3;                                                                 % Maximum effect size displayed on the colorbar
elseif nargin==4
    thres=varargin{1};
    maxd=3;
elseif nargin==5
    thres=varargin{1};
    maxd=varargin{2};
end

if exist('cms_d.mat','file')==2
    load cms_d.mat cms_d
    cms_d=brighten(cms_d,.4);
    if maxd~=3
        cms_d=interp1(1:300,cms_d,linspace(1,300,maxd*100),'spline');
    end
else
    cms_d=colormap(turbo(maxd*100));
    cms_d=brighten(cms_d,.8);
end

for c=1:SPMi.nClusters

    b_inf=round(SPMi.clusters{1,c}.endpoints(1,1));
    b_sup=round(SPMi.clusters{1,c}.endpoints(1,2));
    if b_sup-b_inf < thres
        continue
    end
    d=zeros(1,b_sup-b_inf+1);
    ind=1;

    for b=b_inf+1:b_sup+1

        mY0=mean(Y0(:,b));                                                  % Mean of the distribution
        sY0=std(Y0(:,b));                                                   % Standard deviation
        nY0=size(Y0(:,b),1);                                                % Sample size

        mY1=mean(Y1(:,b));
        sY1=std(Y1(:,b));
        nY1=size(Y1(:,b),1);

        s_pooled=sqrt(((sY0^2*(nY0-1))+(sY1^2*(nY1- 1)))/(nY0+nY1-2));
        d(ind)=abs(mY1-mY0)/s_pooled;                                       % Effect size (Cohen's d)
        ind=ind+1;
    end

    rectangle('Position',[b_inf 0 b_sup-b_inf 180],'EdgeColor','none','FaceColor',cms_d(round(mean(d)*100),:))

end

end