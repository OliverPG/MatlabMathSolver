%·ÉÉ³×ßÊ¯·¨(FeiShaZouShiFa)
%Random Running Stones
%
clear;clc;tic;
clf;
% close all;
syms x y zv
r2=x.^2+2.*y.^2;
z=sqrt(r2)+30;
% z=10-10.*exp(-r2./100).*cos(2.*pi./5.*sqrt(r2));
zFunc=matlabFunction(z);
fFunc=matlabFunction(zFunc-zv);varsList=symvar(sym(fFunc));
xv=linspace(-20,20,100);yv=xv;
[xGrid,yGrid]=ndgrid(xv,yv);
zGrid=zFunc(xGrid,yGrid);
rpics=2;cpics=2;
figure(1);
subplot(rpics,cpics,1);
plot3D(xGrid,yGrid,zGrid,zGrid);xlabel('x');ylabel('y');zlabel('z');hold on;
dt=0.05;%time evolution factor
m0=1;
v0=zeros(1,3);vError=1e-10;err=1e-20;
dzError=1e-3;fError=1e-6;
vmVec=v0;
vmVnorm=vecnorm(vmVec,2,2);
g=9.8;
mu=0.1;%friction factor
nDim=3;
GVec=zeros(1,nDim);GVec(end)=-1;
GVec=GVec.*m0*g;
sVec=jacobian(fFunc,varsList);
sBV=sVec./norm(sVec);
sBVFunc=matlabFunction(sBV);
boundaryLim=...
    [-20,-20,-20;%Lower Lim
      20,20,20];%Higher Lim
  %Dim1,Dim2,Dim3

nPoints=1000;pos=rand(nPoints,nDim);%each row as a point
Ndz=20;dzMat=ones(nPoints,Ndz).*10;
stopFlag=0;activePList=1:nPoints;
pos=repmat(boundaryLim(1,:),nPoints,1)+repmat(boundaryLim(2,:)-boundaryLim(1,:),nPoints,1).*pos;
pos(activePList,end)=zFunc(pos(activePList,1),pos(activePList,2));%Initial Positions
GVec=repmat(GVec,nPoints,1);
vmVec=repmat(vmVec,nPoints,1);
vmVnorm=vecnorm(vmVec,2,2);vmVnorm=repmat(vmVnorm,1,nDim)+err;
dtVec=ones(nPoints,1).*dt;
count=0;
while stopFlag==0 && count<400 || count==0
    posActiCell=mat2cell(pos(activePList,1:nDim-1),length(activePList),ones(1,nDim-1));
    sBV_v=sBVFunc(posActiCell{:});
%     sBV_v=sBVFunc(pos(activePList,1),pos(activePList,2));
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
    Ev=m0.*vmVnorm(:,1).^2./2;
    rVec(activePList,:)=vmVec(activePList,:).*dtVec(activePList,:);
    oldPos=pos;
    pos(activePList,:)=pos(activePList,:)+rVec(activePList,:);
    pos(activePList,end)=zFunc(pos(activePList,1),pos(activePList,2));
    dz=pos(:,end)-oldPos(:,end);
    dzMat=circshift(dzMat,-1,2);dzMat(:,end)=dz;dzMatAver=mean(dzMat,2);
    figure(1);subplot(rpics,cpics,1);
    if exist('handleS','var')
        delete(handleS);
    end
    hold on;
    handleS=scatter3(pos(:,1),pos(:,2),pos(:,3),'.');title('POS');
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
    [owR,owC]=find(pos(:,1:end-1)-repmat(boundaryLim(1,1:end-1),nPoints,1)<0);
    [perR,perC]=find(pos(:,1:end-1)-repmat(boundaryLim(2,1:end-1),nPoints,1)>0);
    overIndex=unique([owR;perR]);
    if ~isempty(overIndex)
        pos(overIndex,:)=oldPos(overIndex,:);
    end
%     stopIndex=unique([v0Index;dz0Index;dzA0Index;overIndex]);
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
scatter3(pos(:,1),pos(:,2),pos(:,3),'.')
% axis equal
% vmV=feval(vmVec,pos(:,1),pos(:,2))
toc;