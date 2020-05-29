%%%%%%%%%%%%%%%%%%%%%%%%%
%y(x)=Gaussian Funtion
%y(0)=0,y'(0)=k1
%y(L)=-h0,y'(L)=k2
%%%%%%%%%%%%%%%%%%%%%%%%%
%   4 for
%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
output=1;%output txt file or not
tic
syms y x h b c d
yi(x)=exp(x.^2);yip(x)=diff(yi,x);
y=h.*exp(-((x-b)./c).^2./2)+d;
syms L h0 k1 k2
syms f1
count=0;countSkip=0;skipParaMat=[];
f1=h0+d.*(1-exp(-L.*(k2+k1.*(1+h0./d))./(2.*(h0+d))));
f1A=1+h0./d;
f1B=exp(-L.*(k2+k1.*(1+h0./d))./(2.*(h0+d)));
%%%Parameters set%%%%%%%%%%%%
x0=0;y0=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0_v=500./2;
L_v=linspace(150,1000,5);
L_v=(L_v(1:end-1)+L_v(2:end))./2;
angle1=0; k1_v=tand(-angle1);
angle2=linspace(1,45,5);
angle2=(angle2(1:end-1)+angle2(2:end))./2;
k2_v=tand(-angle2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalCount=length(h0_v)*length(L_v)*length(angle1)*length(angle2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
for h0_vi=h0_v
    for L_vi=L_v
        for k1_vi=k1_v
            for k2_vi=k2_v
                fprintf('\n######################################################');
                fprintf('\nParameters:\n\th0=%.4f\n\tL=%.4f\n\tk1=%.4f\n\tk2=%.4f',h0_vi,L_vi,k1_vi,k2_vi);
                f1_d(d)=vpa(subs(f1,[L,h0,k1,k2],[L_vi,h0_vi,k1_vi,k2_vi]));
                D2=-h0_vi;D1=-k1_vi.*h0_vi./(k1_vi+k2_vi);
                switch k1_vi+k2_vi>0
                    case true
                        d_start=(D1+D2)./2;
                    case false
                        d_start=[1.5.*D2-eps(D2),D1./2+eps(D1),eps(1)];
                end
                d_v=vpasolve(f1_d,d_start(1));
                deltaD=vpa(f1_d(d_v));
                if isempty(deltaD)
                    deltaD=1;
                end
                countLoop=1;
                while ~(abs(deltaD)<(1e-10) || countLoop==3)
                    countLoop=countLoop+1;
                    d_v=vpasolve(f1_d,d_start(countLoop));
                    deltaD=vpa(f1_d(d_v));
                    if isempty(deltaD)
                        deltaD=1;
                    end
                    fprintf('\ncountLoop=%d\tdeltaD=%f',countLoop,deltaD);
                end
                tryLeft=D2-1000;tryRight=D2+1000;
                definput = {'-500','0.1'};%countLoop=3;                
                f1A_d(d)=vpa(subs(f1A,h0,h0_vi));
                f1B_d(d)=vpa(subs(f1B,[L,h0,k1,k2],[L_vi,h0_vi,k1_vi,k2_vi]));
                while abs(deltaD)>(1e-10)
                    figure(2);subplot(2,2,1);
                    plotFuncOn(tryLeft,tryRight,f1B_d,2);
                    subplot(2,2,2);
                    plotFuncOn(tryLeft,tryRight,f1A_d,2);
                    subplot(2,2,[3,4]);
                    plotFuncOn(tryLeft,tryRight,f1B_d,2);
                    hold on;plotFuncOn(tryLeft,tryRight,f1A_d,2);
                    set(gcf,'outerposition',get(0,'screensize'));
                    prompt = {'x Start:','x End:'};
                    title = 'Input Scale ';
                    dims = [1 35];
                    xCell = inputdlg(prompt,title,dims,definput);
                    tryLeft=str2double(xCell{1});tryRight=str2double(xCell{2});
                    definput = {num2str(tryLeft+(tryRight-tryLeft)./10),num2str(tryRight-(tryRight-tryLeft)./10)};
                    d_v=vpasolve(f1_d,[tryLeft,tryRight]);
                    deltaD=vpa(f1_d(d_v));
                    if isempty(deltaD)
                        deltaD=1;
                    else                        
                    end
                    countLoop=countLoop+1;
                    fprintf('\ncountLoop=%d\tdeltaD=%f',countLoop,deltaD);
                    close(2);
%                     close(3);
                end
                if abs(deltaD)<(1e-10)
                    c_v=sqrt(L_vi.*(h0_vi+d_v)./(k2_vi-k1_vi.*(1+h0_vi./d_v)));
                    h_v=-d_v.*exp(c_v.^2.*k1_vi.^2./(2.*d_v.^2));
                    b_v=-c_v.^2.*k1_vi./d_v;
                    ySolver=vpa(subs(y,[h,b,c,d],[h_v,b_v,c_v,d_v]));
                    fprintf('\nResult:\n\th=%.4f\n\tb=%.4f\n\tc=%.4f\n\td=%.4f',h_v,b_v,c_v,d_v);
                    fprintf('\ny=%s',ySolver);
                    yS(x)=ySolver;yp(x)=diff(yS,x);
                    x1=x0+L_vi;y1=y0-h0_vi;
                    RCheck=abs([double(yS(x0)-y0),double(yS(x1)-y1),double(yp(x0)-k1_vi),double(yp(x1)-k2_vi)]);
                    if sum(RCheck<1e-10)~=4
                        fprintf('\nResult Check Wrong:\n\tyS(x0)-y0=%f\n\tyS(x1)-y1=%f\n\typ(x0)-k1=%f\n\typ(x1)-k2=%f',double(yS(x0)-y0),double(yS(x1)-y1),double(yp(x0)-k1_vi),double(yp(x1)-k2_vi));
                        skipParaMat=[skipParaMat;L_vi,h0_vi,k1_vi,k2_vi];
                        continue;
                    end
                    %                     fprintf('\nResult Check:\n\tyS(x0)-y0=%f\n\tyS(x1)-y1=%f\n\typ(x0)-k1=%f\n\typ(x1)-k2=%f',double(yS(x0)-y0),double(yS(x1)-y1),double(yp(x0)-k1_vi),double(yp(x1)-k2_vi));
                    hold on;plotFuncOn(0,L_vi,ySolver,1);
                    hold on;plot([0,L_vi],[0,-h0_vi],'o');
                    axis equal
                    count=count+1;
                    if output
                        paraStr=sprintf('(L%.1fH%.1fA1%.1fA2%.1f)',L_vi,h0_vi,atand(k1_vi),atand(k2_vi));
                        syms t
                        chF1=char(t);%X or U
                        chF2=num2str(0);%Y or V
                        chF3=char(yS(t));%Z or N
                        if count==1
                            dateStr=datestr(now,'yyyymmddHHMMSS');
                            fileName=sprintf('GaussCurve%s.txt',dateStr);
                            fileID=fopen(fileName,'w');
                        end
                        fprintf(fileID,'%s\n%s\n%s\n%s\n',paraStr,chF1,chF2,chF3);
                        fprintf('\nWrite to %s :\n\t%s\n\t%s\n\t%s\n\t%s',fileName,paraStr,chF1,chF2,chF3);
                    end
                    fprintf('\n%d/%d (count/total)',count,totalCount);
                end
            end
        end
    end
end
fclose all;
toc