%·ÉÉ³×ßÊ¯·¨(FeiShaZouShiFa)
%Random Running Stones
%Expulso

clear;clc;
syms a1 a2 a3 a4
% Test part1, formulas for different dims
% r2=x1.^2+2.*x2.^2;
% r3=x1.^2+3.*x2.^2+2.*x3.^2+2.*x4.^2;
% z=cos(x1).*x1;
% z=x1.^5./10000;
% z=sqrt(r2)+30;
% z=10-10.*exp(-r2./100).*cos(2.*pi./5.*sqrt(r2));
% z=sqrt(r3)+30;
% boundaryLim=...
%     [-20,-20,-20,-20,-20;%Lower Lim
%     20,20,20,20,20];%Higher Lim
%Dim1,Dim2,Dim3
% Test part2, formulas for fitting scatters
boundaryLim=...
    [0,-2,-20,-20,-20;%Lower Lim
    5,20,20,20,20];%Higher Lim
syms x yv
yInit=cos(x./pi*4)*10+1;yInitFunc=matlabFunction(yInit);
yFit=cos(x.*a1)*a2+a3;yFitFunc=matlabFunction(yFit);
yf=yFit-yv;
yfFunc=matlabFunction(yf);
nScatters=50;x1Vec=linspace(boundaryLim(1,1),boundaryLim(2,1),nScatters)';
scatters=[x1Vec,yInitFunc(x1Vec)+randi(5,nScatters,1)./5-0.5];
scatterFig=figure;
plot(scatters(:,1),scatters(:,2),'*');
iVars=[4,5];
yfFuncs=funcByVal(yfFunc,iVars,scatters);
% varsV=num2cell(symvar(yf));varsV(iVars)=mat2cell(scatters,length(scatters),[1,1]);
% z=sqrt(simplify(sum((yfFuncs./sqrt(1+x1.^2)).^2))./(1e1));
% z=sqrt(simplify(sum((yfFuncs./1).^2))./(1e1));
z=sqrt(sum((yfFuncs./1).^2)./(1e1));
zFunc=matlabFunction(z);
minPos=localMinsRRS1(zFunc,boundaryLim);
[minZ,minR]=min(minPos(:,end));
yR=funcByVal(yFitFunc,[1:3],minPos(minR,1:3));
% hold on;
plotFuncOn(0,5,matlabFunction(yR),scatterFig);

function [pos,fig]=localMinsRRS1(zFunc,boundaryLim)
tic;
syms zv
fFunc=matlabFunction(zFunc-zv);varsList=symvar(sym(fFunc));varsListStr=string(varsList);
nVar=length(varsList);

rpics=2;cpics=2;
fig=figure;
subplot(rpics,cpics,1);
if nVar>=3
    x1v=linspace(boundaryLim(1,1),boundaryLim(2,1),100);
    x2v=linspace(boundaryLim(1,2),boundaryLim(2,2),100);
    [x1Grid,x2Grid]=ndgrid(x1v,x2v);
    xyuGridCell=[{x1Grid},{x2Grid}];
    if nVar>3
        uCell=num2cell(mean(boundaryLim(:,3:nVar-1)));
        xyuGridCell=[xyuGridCell,uCell];
    end
    zGrid=zFunc(xyuGridCell{:});
    plot3D(x1Grid,x2Grid,zGrid,zGrid);xlabel('Dim1');ylabel('Dim2');zlabel('z');hold on;
else % nVar==2
    zGrid=zFunc(x1v);
    plot(x1v,zGrid,'--','LineWidth',0.5,'color','black');hold on;
end
err=1e-20;% zero
vError=1e-10;% zero for velocity
dzError=1e-3;% zero for dz
nDim=nVar;% motion dimision or number of variates
dt=0.05;% time evolution factor
m0=1;
g=9.8;
mu=0.2;%friction factor
GVec=zeros(1,nDim);GVec(end)=-1;
GVec=GVec.*m0*g;
sVec=jacobian(fFunc,varsList);
sBV=sVec./norm(sVec);
sBVFunc=matlabFunction(sBV);varBvList=symvar(sBV);varBvListStr=string(varBvList);
[BvIndex,varsIndex]=find(varsListStr==varBvListStr');nVarBv=length(varBvList);

nPoints=100;
activePList=1:nPoints;
pos=rand(nPoints,nDim);% Each row as the position of a point
pos=repmat(boundaryLim(1,1:nDim),nPoints,1)+repmat(boundaryLim(2,1:nDim)-boundaryLim(1,1:nDim),nPoints,1).*pos;
posActiCell_z=mat2cell(pos(activePList,1:nDim-1),length(activePList),ones(1,nDim-1));
pos(activePList,end)=zFunc(posActiCell_z{:});% Initial positions on surface
GVec=repmat(GVec,nPoints,1);
vmVec=zeros(nPoints,nDim);
vmVnorm=zeros(nPoints,nDim)+err;
dtVec=ones(nPoints,1).*dt;
Ndz=20;dzMat=ones(nPoints,Ndz).*10;
stopFlag=0;
count=0;
countA=0;
while stopFlag==0 && count<400 || countA<=Ndz
    if 0
    % Points fission
    rFiss=min(boundaryLim(2,:)-boundaryLim(1,:))/(nPoints/1);
    nFiss=9;%Each Dim
    [pFissList,pFissList3D]=pointsOnCubeFace(pos(activePList,1:end-1),rFiss,nFiss);
    nTotalFiss=length(pFissList(:,1));
    [r_pFiss,c_pFiss,p_pFiss]=size(pFissList3D);
    nTotalFiss=r_pFiss*p_pFiss;
    % Check over boundary
    pFissActiveList=1:nTotalFiss;
    expandLowBd=repmat(boundaryLim(1,1:nDim-1),r_pFiss,1,p_pFiss);
    expandUpBd=repmat(boundaryLim(2,1:nDim-1),r_pFiss,1,p_pFiss);
    indexLow=find(pFissList3D-expandLowBd<0);
    indexUp=find(pFissList3D-expandUpBd>0);
    if ~isempty(indexLow)
        pFissList3D(indexLow)=expandLowBd(indexLow);
    end
    if ~isempty(indexUp)
        pFissList3D(indexUp)=expandUpBd(indexUp);
    end
    pFissCell=mat2cell(pFissList3D,r_pFiss,ones(1,c_pFiss),p_pFiss);
    zFissList=reshape(zFunc(pFissCell{:}),r_pFiss,p_pFiss);
    [zFissMinList,zFissMinIndex]=min(zFissList,[],2);
    indexMove=find(zFissMinList<pos(activePList,end));
    indexMove=[];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(indexMove)
        nMove=length(indexMove);
        select_pFissList=[];
        for iMove=1:nMove
            select_pFissList(iMove,:)=pFissList3D(indexMove(iMove),:,zFissMinIndex(iMove));
        end
        %         copy_pFissList=pFissList3D(indexMove,:,zFissMinIndex);
        oldPos=pos;
        pos(activePList(indexMove),:)=[select_pFissList,zFissMinList(indexMove)];
        vmVec(activePList(indexMove),:)=zeros(nMove,nDim);
        count=count+1;
    end
    end
    
    indexMove=[];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Points falling
    [activePList_noMove,activePList_noMoveIndex]=setdiff(activePList,activePList(indexMove));
    if 1
    posActiCell_BV=mat2cell(pos(activePList,varsIndex),length(activePList),ones(1,nVarBv));
    sBV_v=sBVFunc(posActiCell_BV{:});
    FgpVec(activePList,:)=GVec(activePList,:)-repmat(dot(GVec(activePList,:),sBV_v,2),1,nDim).*sBV_v;
    fVec=zeros(size(FgpVec));
%     if ~isempty(activePList_noMove)
%         fVec(activePList_noMove,:)=-abs(repmat(dot(GVec(activePList_noMove,:),sBV_v(activePList_noMoveIndex,:),2),1,nDim)).*vmVec(activePList_noMove,:)./vmVnorm(activePList_noMove,:).*mu;
%     end
    fVec(activePList,:)=-abs(repmat(dot(GVec(activePList,:),sBV_v(:,:),2),1,nDim)).*vmVec(activePList,:)./vmVnorm(activePList,:).*mu;
    aVec(activePList,:)=(FgpVec(activePList,:)+fVec(activePList,:))./m0;
    dvVec(activePList,:)=aVec(activePList,:).*dtVec(activePList,:);
    vmVec(activePList,:)=vmVec(activePList,:)+dvVec(activePList,:);
    vmVnorm=vecnorm(vmVec,2,2);vmVnorm=repmat(vmVnorm,1,nDim)+err;
    rVec(activePList,:)=vmVec(activePList,:).*dtVec(activePList,:);
    oldPos=pos;
    pos(activePList,:)=pos(activePList,:)+rVec(activePList,:);
    posActiCell_z=mat2cell(pos(activePList,1:nDim-1),length(activePList),ones(1,nDim-1));
    pos(activePList,end)=zFunc(posActiCell_z{:});
    countA=count+1;
    end
    
    dz=pos(:,end)-oldPos(:,end);
    dzMat=circshift(dzMat,-1,2);dzMat(:,end)=dz;dzMatAver=mean(dzMat,2);
    figure(fig);subplot(rpics,cpics,1);
    if exist('handleS','var')
        delete(handleS);
    end
    hold on;
    if nVar>=3
        handleS=scatter3(pos(:,1),pos(:,2),pos(:,end),'.');title('Position');
    else %nvar=2
        handleS=scatter(pos(:,1),pos(:,2),'.');title('Position');
    end
    subplot(rpics,cpics,2);hold on;scatter(repmat(count,nPoints,1),pos(:,end),'.');title('z');
    if count>Ndz
        subplot(rpics,cpics,3);hold on;scatter(repmat(count,nPoints,1),dzMatAver,'.');title('dzA');
    end
    subplot(rpics,cpics,4);hold on;scatter(repmat(count,length(activePList),1),activePList,'.');title('activePList');
    drawnow;
    %velocity near zero
    v0Index=find(abs(vmVec(:,end))<vError);
    %dz near zero
    dz0Index=find(abs(dz)<dzError);
    %dzNAver near zero
    dzA0Index=find(abs(dzMatAver)<dzError);
    %point over boundary
    [owR,owC]=find(pos(:,1:end-1)-repmat(boundaryLim(1,1:nDim-1),nPoints,1)<0);
    [perR,perC]=find(pos(:,1:end-1)-repmat(boundaryLim(2,1:nDim-1),nPoints,1)>0);
    overIndex=unique([owR;perR]);
    if ~isempty(overIndex)
        pos(overIndex,:)=oldPos(overIndex,:);
    end
    stopIndex=unique([dzA0Index;overIndex]);
%     stopIndex=unique([v0Index;dzA0Index;overIndex]);
    if ~isempty(stopIndex)
        activePList=setdiff(activePList,stopIndex);
    end
    if isempty(activePList)
        stopFlag=1;
    end
    %Points Rising control
    riseIndex=find(dz>0);
    if ~isempty(riseIndex)
        vmVec(riseIndex,:)=zeros(length(riseIndex),nDim);
        vmVnorm=vecnorm(vmVec,2,2);vmVnorm=repmat(vmVnorm,1,nDim)+err;
        %         pos=oldPos;
    end
    count=count+1;
%     end
end
toc;
end