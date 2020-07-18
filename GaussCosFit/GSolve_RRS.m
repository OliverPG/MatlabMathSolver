%��ɳ��ʯ��(FeiShaZouShiFa)
%Random Running Stones
%Expulso

clear;clc;tic;
clf;
syms x1 x2 x3 x4
r2=x1.^2+2.*x2.^2;
r3=x1.^2+3.*x2.^2+2.*x3.^2+2.*x4.^2;
% z=cos(x1).*x1;
% z=x1.^5./10000;
% z=sqrt(r2)+30;
% z=10-10.*exp(-r2./100).*cos(2.*pi./5.*sqrt(r2));
z=sqrt(r3)+30;
zFunc=matlabFunction(z);
boundaryLim=...
    [-20,-20,-20,-20,-20;%Lower Lim
      20,20,20,20,20];%Higher Lim
  %Dim1,Dim2,Dim3
syms zv
fFunc=matlabFunction(zFunc-zv);varsList=symvar(sym(fFunc));varsListStr=string(varsList);
nVar=length(varsList);

% xv=linspace(num2cell(:,1),100);
x1v=linspace(boundaryLim(1,1),boundaryLim(2,1),100);
rpics=2;cpics=2;
figure(1);
subplot(rpics,cpics,1);
if nVar>=3
    x2v=linspace(boundaryLim(1,2),boundaryLim(2,2),100);
    [x1Grid,x2Grid]=ndgrid(x1v,x2v);
    xyuGridCell=[{x1Grid},{x2Grid}];
    if nVar>3
        uCell=num2cell(mean(boundaryLim(:,3:end-1)));
        xyuGridCell=[xyuGridCell,uCell];
    end
    zGrid=zFunc(xyuGridCell{:});
    plot3D(x1Grid,x2Grid,zGrid,zGrid);xlabel('x');ylabel('y');zlabel('z');hold on;
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

nPoints=1000;
activePList=1:nPoints;
pos=rand(nPoints,nDim);% Each row as the position of a point
pos=repmat(boundaryLim(1,nDim),nPoints,1)+repmat(boundaryLim(2,nDim)-boundaryLim(1,nDim),nPoints,1).*pos;
posActiCell_z=mat2cell(pos(activePList,1:nDim-1),length(activePList),ones(1,nDim-1));
pos(activePList,end)=zFunc(posActiCell_z{:});% Initial positions on surface
GVec=repmat(GVec,nPoints,1);
vmVec=zeros(nPoints,nDim);
vmVnorm=zeros(nPoints,nDim)+err;
dtVec=ones(nPoints,1).*dt;
Ndz=20;dzMat=ones(nPoints,Ndz).*10;
stopFlag=0;
count=0;
while stopFlag==0 && count<400 || count==0
    posActiCell_BV=mat2cell(pos(activePList,varsIndex),length(activePList),ones(1,nVarBv));
    sBV_v=sBVFunc(posActiCell_BV{:});
    FgpVec(activePList,:)=GVec(activePList,:)-repmat(dot(GVec(activePList,:),sBV_v,2),1,nDim).*sBV_v;
    if count==0
        fVec=zeros(size(FgpVec));
        fOld=fVec;
    else
        fOld=fVec;
        fVec(activePList,:)=-abs(repmat(dot(GVec(activePList,:),sBV_v,2),1,nDim)).*vmVec(activePList,:)./vmVnorm(activePList,:).*mu;
    end
    aVec=(FgpVec+fVec)./m0;
    dvVec=aVec.*dtVec;
    vmVec(activePList,:)=vmVec(activePList,:)+dvVec(activePList,:);
    vmVnorm=vecnorm(vmVec,2,2);vmVnorm=repmat(vmVnorm,1,nDim)+err;
    rVec(activePList,:)=vmVec(activePList,:).*dtVec(activePList,:);
    oldPos=pos;
    pos(activePList,:)=pos(activePList,:)+rVec(activePList,:);
    posActiCell_z=mat2cell(pos(activePList,1:nDim-1),length(activePList),ones(1,nDim-1));
    pos(activePList,end)=zFunc(posActiCell_z{:});
    dz=pos(:,end)-oldPos(:,end);
    dzMat=circshift(dzMat,-1,2);dzMat(:,end)=dz;dzMatAver=mean(dzMat,2);
    figure(1);subplot(rpics,cpics,1);
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
    stopIndex=unique([v0Index;dzA0Index;overIndex]);
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
end
% minPos=pos(find(abs(pos(:,end)-min(pos(:,end)))<0.1),:);
% figure;plot(minPos(:,1),minPos(:,2),'*');toc;
figure;
% plot(pos(:,1),pos(:,2),'*');
if nVar>=3
    scatter3(pos(:,1),pos(:,2),pos(:,end),'.');
else % nVar=2
    scatter(pos(:,1),pos(:,2),'.');
end
% axis equal
% vmV=feval(vmVec,pos(:,1),pos(:,2))
toc;