%%%%%%%%%%%%%%%%%%%%%%%%%
%y(x)=Gaussian Funtion
%y(0)=0,y'(0)=k1
%y(L)=-h0,y'(L)=k2
%%%%%%%%%%%%%%%%%%%%%%%%%
%   4 for
%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
format long
err=1e-9;
InfN=1e8;
output=1;%output txt file or not
tic
syms y x h b c d
yi(x)=exp(x.^2);yip(x)=diff(yi,x);
y=h.*exp(-((x-b)./c).^2./2)+d;
syms L h0 k1 k2
syms f1 f1_d
count=0;unsolveCount=0;unsolvedPara=[];
f1=h0+d.*(1-exp(-L.*(k2+k1.*(1+h0./d))./(2.*(h0+d))));
f1A=1+h0./d;
f1B=exp(-L.*(k2+k1.*(1+h0./d))./(2.*(h0+d)));
%%%Parameters set%%%%%%%%%%%%
x0=0;y0=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0_v=500./2;
L_v=linspace(150,1000,9);
% L_v=(L_v(1:end-1)+L_v(2:end))./2;
angle1=0; k1_v=tand(-angle1);
angle2=linspace(1,45,9);
% angle2=(angle2(1:end-1)+angle2(2:end))./2;
k2_v=tand(-angle2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalCount=length(h0_v)*length(L_v)*length(angle1)*length(angle2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figHd=figure;
for h0_vi=h0_v
    for L_vi=L_v
        for k1_vi=k1_v
            for k2_vi=k2_v
                fprintf('\n######################################################');
                fprintf('\nParameters:\n\th0=%.4f\n\tL=%.4f\n\tk1=%.4f\n\tk2=%.4f',h0_vi,L_vi,k1_vi,k2_vi);
                vars=[L_vi,h0_vi,k1_vi,k2_vi];
                vars=sym(vars*eps)/eps;
                f1_d(d)=vpa(subs(f1,[L,h0,k1,k2],vars));
                D2=-h0_vi;D1=double(-k1_vi.*h0_vi./(k1_vi+k2_vi));
                d0Vec=[];
                if k1_v+k2_v<0
                    [d0Vec(1,1),d0Vec(1,2)]=solveMono(f1_d,[-InfN,D2]);
                    if abs(d0Vec(1,2))>err || abs(d0Vec(1,1)+h0_vi)<err
                        [d0Vec(2,1),d0Vec(2,2)]=solveMono(f1_d,[0,InfN]);
                        if abs(d0Vec(2,2))>err && D1~=0
                            [d0Vec(3,1),d0Vec(3,2)]=solveMono(f1_d,[D1,0]);
                        end
                    end
                elseif k1_v+k2_v>0
                    [d0Vec(1,1),d0Vec(1,2)]=solveMono(f1_d,[D2,D1]);
                else
                    d0Vec(1,:)=[D2,f1_d(D2)];
                end
                d_v=d0Vec(end,1);
                deltaD=d0Vec(end,2);
                if abs(deltaD)<err
                    c_v=sqrt(L_vi.*(h0_vi+d_v)./(k2_vi-k1_vi.*(1+h0_vi./d_v)));
                    h_v=-d_v.*exp(c_v.^2.*k1_vi.^2./(2.*d_v.^2));
                    b_v=-c_v.^2.*k1_vi./d_v;
                    solVars=[h_v,b_v,c_v,d_v];
                    solVars=sym(solVars*eps)/eps;
                    ySolver=vpa(subs(y,[h,b,c,d],solVars));
                    fprintf('\nResult:\n\th=%.4f\n\tb=%.4f\n\tc=%.4f\n\td=%.4f',h_v,b_v,c_v,d_v);
                    fprintf('\ny=%s',ySolver);
                    yS(x)=ySolver;yp(x)=diff(yS,x);
                    x1=x0+L_vi;y1=y0-h0_vi;
                    RCheck=abs([double(yS(x0)-y0),double(yS(x1)-y1),double(yp(x0)-k1_vi),double(yp(x1)-k2_vi)]);
                    if sum(RCheck<err)~=4
                        fprintf('\nResult Check Wrong:\n\tyS(x0)-y0=%f\n\tyS(x1)-y1=%f\n\typ(x0)-k1=%f\n\typ(x1)-k2=%f',double(yS(x0)-y0),double(yS(x1)-y1),double(yp(x0)-k1_vi),double(yp(x1)-k2_vi));
                        unsolveCount=unsolveCount+1;
                        unsolvedPara(unsolveCount,:)=[L_vi,h0_vi,k1_vi,k2_vi];
                        continue;
                    end
                    hold on;plotFuncOn(0,L_vi,ySolver,figHd);
                    hold on;plot([0,L_vi],[0,-h0_vi],'o');
                    axis equal
                    count=count+1;
                    if output
                        paraStr=sprintf('(L%gH%gA1%gA2%g)',L_vi,h0_vi,atand(k1_vi),atand(k2_vi));
                        syms t
                        chF1=char(t);%X or U
                        chF2=num2str(0);%Y or V
                        chF3=char(yS(t));%Z or N
                        if count==1
                            dateStr=datestr(now,'yyyymmdd_HHMMSS');
                            fileName=sprintf('GaussCurve(%d)%s.txt',totalCount,dateStr);
                            fileID=fopen(fileName,'w');
                        end
                        fprintf(fileID,'%s\n%s\n%s\n%s\n',paraStr,chF1,chF2,chF3);
                        fprintf('\nWrite to %s :\n\t%s\n\t%s\n\t%s\n\t%s',fileName,paraStr,chF1,chF2,chF3);
                    end
                    fprintf('\n%d/%d (count/total)',count,totalCount);
                else
                    unsolveCount=unsolveCount+1;
                    unsolvedPara(unsolveCount,:)=[L_vi,h0_vi,k1_vi,k2_vi];
                end
            end
        end
    end
end
fclose all;
toc