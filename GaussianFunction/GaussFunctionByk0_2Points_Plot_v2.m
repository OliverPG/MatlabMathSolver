%1-D gauss distribution 
%由两个条件确定一条高斯曲线：
%1 原点处的斜率为k0 （过原点）
%2 左半支曲线过(Lx,hy)
%画出这个线段
clear;clc;
tic
error=0;
errorI=0;errorStore=0;
count=0;

theta=5.6014114;%deg half angle of the vertex
% LCut=500;
wCut=106.181451*2;%LCut*tand(theta)*2; 
wCutHalf=wCut/2;
k0=tand(theta);
k1_mat=[1,0.9,0.1];%[1,0.98,0.9,0.8,0.6,0.5,0.3,0.1,0.01];%,0.01,1/2 %大于0 小于等于1 为1时（Lx,hy）为曲线顶点
Lx_mat=142;%linspace(LCut/2,wCutHalf./2,3);%线段的终点x坐标
hy=wCutHalf;%线段的终点y坐标
[~,cs1]=size(k1_mat);
[~,cs2]=size(Lx_mat);
nColors=colormap(jet(cs1*cs2));
syms vh vMean vSigm vd
X=[vh,vMean,vSigm,vd];
nPoints=1000;
for Lx=Lx_mat
    for k1=k1_mat
        count=count+1
        vh0=wCut;
        vMean0=Lx;%./(1-sqrt((h-hy)./(2.*hy.*log(2))))        
        vSigm0=vMean0./sqrt(2.*log(2));
        vd0=0;
        x0=[vh0,vMean0,vSigm0,vd0];
        func=@(X)formulasV2(X,k0,k1,hy,Lx);
        tolerance=1e-12;
        options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FunctionTolerance',tolerance,'OptimalityTolerance',tolerance,'StepTolerance',tolerance,'MaxIterations',800,'PlotFcn',@optimplotfirstorderopt);
%         options = optimoptions('fsolve','Display','none','Algorithm','trust-region-reflective','PlotFcn',@optimplotfirstorderopt);
        %'levenberg-marquardt' 'trust-region-dogleg' 'trust-region'
        [X_solver,fval,exitflag,output]=fsolve(func,x0,options);output
        if norm(fval)>1e-3
            errorI=errorI+1;
            errorStore(errorI,1)=count;
            errorStore(errorI,2:5)=fval;
        end
        h=X_solver(1);mean=X_solver(2);sigm=X_solver(3);d=X_solver(4);
        
        x=linspace(0,Lx,nPoints);
        y=h.*exp(-(x-mean).^2./(2.*sigm.^2))+d;  
        y0_d(count,1)=-(0-mean)./(sigm.^2).*h.*exp(-(0-mean).^2./(2.*sigm.^2));
        figure(1)
        plot(x,y,'color',nColors(count,:),'linewidth',2);
        hold on; 
        drawnow;
        legendStr{count}=sprintf('%d k1=%f,Lx=%f',count,k1,Lx);
        syms vx
%         vy=vpa(h*exp(-(vx-mean)^2/(2*sigm^2))+d);  
        vy=(h*exp(-(vx-mean)^2/(2*sigm^2))+d);  
        curData{count,1}=sprintf('k1=%f\tLx=%f\thy=%f',k1,Lx,hy);
        curData{count,2}=sprintf('xStart\txEnd\tyStart\tyEnd');
        curData{count,3}=sprintf('%f\t%f\t%f\t%f',x(1),x(end),y(1),y(end));
        curData{count,4}=char(vy);
        curData{count,5}=sprintf('__________________________________%d',count);
        
%         if Lx~=Lx_mat(2)
%             break;
%         end
    end
end
if error==0
    fileID=fopen('GaussFunctionPlot_Cache.txt','w');
    for jj=1:1:count
        for ii=1:1:5
            fprintf(fileID,'%s\n',curData{jj,ii});
        end
    end
    timeStr=datetime;
    fprintf(fileID,'%s\n',timeStr);    
    fclose(fileID);
    type GaussFunctionPlot_Cache.txt
    hold on;
%     xWedg=linspace(0,wCutHalf,1000);
    xWedg=linspace(0,Lx,1000);
    kWedg=wCutHalf./Lx;
    yWedg=kWedg.*xWedg;
    plot(xWedg,yWedg,'--','color','black','linewidth',2);
    count=count+1;
    legendStr{count}=sprintf('%d k_W_e_d_g_e=%f',count,kWedg);
    legend(legendStr,'Location','bestoutside');
    xlabel('x(mm)')
    ylabel('y(mm)')
    title('Gaussian Function');
    axis([0 max(Lx_mat) 0 hy]);
    grid on
    errorStore
    atand(y0_d)
    dispStr=sprintf('Done');
end
toc