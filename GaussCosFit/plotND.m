clear;clc;clf;
nPoints=10;
nDimm=9;
posNd=randi(30,nPoints,nDimm);% One point at each row

[rsPos,nDim]=size(posNd);
coordTheta=pi/nDim;
vecBs=[cos((0:nDim-1).*coordTheta)',sin((0:nDim-1).*coordTheta)'];% One base vector at echo row
pos2D=posNd*vecBs;
R=30;
startPos=zeros(nDim,2);
endPos=30.*vecBs;
xPos=[startPos(:,1)';endPos(:,1)'];%[x1;x2]
yPos=[startPos(:,2)';endPos(:,2)'];%[y1;y2]
subplot(2,1,1);
line(xPos,yPos);% Line at each Column
hold on;scatter(pos2D(:,1)',pos2D(:,2)');
axis equal
subplot(2,1,2);
if nDim>=3
    scatter3(posNd(:,1)',posNd(:,2)',posNd(:,3)');
elseif nDim==2
    scatter(posNd(:,1)',posNd(:,2)');
end
    