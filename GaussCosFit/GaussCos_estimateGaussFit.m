clear;clc;
% syms x y 
% syms a1 a2 a3 a4 a5 f00
% a=[a1,a2,a3,a4,a5];
% f=a1*exp(a2*x^2+a3*x)*cos(a4*x)+a5-y;
% f=a1+a2+a3+a4+a5;
% for indexa=1:length(a)
%     b{indexa}=solve(f,a(indexa));
% end
fre=[300;470;640;725;810;980;1150;1320;1490;1575;1660;1830;2000];
rcs=[3.50266273829231;7.18167481560552;8.29853604322728;8.27168802473965;7.23349089331251;4.78699634777701;-0.446225679851643;-10.2257061777637;-6.08043729857107;-3.19274369835503;-1.62180305353850;-0.425360901742028;-1.50920384069627];
a7Esti=mean(rcs);
gaussY=rcs-a7Esti;
gausHd=@(x)sum((x(1).*exp(-x(2).*(fre-x(3)).^2)+x(4)).^2);
% fsHd=@(x)sum((x(1).*exp(-1./x(2).*(fre-x(3)).^2).*cos(x(4).*fre)+x(5)-rcs).^2);
[x0,fval,exitFlag]=GlobalSolve(gausHd,4)
plot(fre,gaussY,'.')
syms xx
ty=x0(1).*exp(-x0(2).*(xx-x0(3)).^2)+x0(4);
%%
% ty=10.*exp(-1/4e5.*(xx+200).^2).*cos(0.5/pi.*(xx+150)./10)+5;%;
hold on;plotFuncOn(300,2000,ty,1);
%%
% x0 = particleswarm(funcHd1,5,-inf(1,5),inf(1,5),options);
% GlobalSolve(funcHd,5)
% GlobalSolve(fs,a,5)
% fsolve(funcHd,[1,1,1,1,1])
% func1(a,b,c,d)=y(p1(1))-p1(2);
% func2(a,b,c,d)=y(p2(1))-p2(2);
% func3(a,b,c,d)=dy(p1(1))-k1;
% func4(a,b,c,d)=dy(p2(1))-k2;
% funcs=[func1(a,b,c,d);func2(a,b,c,d);func3(a,b,c,d);func4(a,b,c,d)];
% vs=[a,b,c,d];
% vs0=[h,L,1.5*L/2,h];
% funcHd=matlabFunction(funcs,'vars',{[a,b,c,d]});
% n=4;m=20;
% [vss,fval,exitFlag] = GlobalSolve(funcHd,5);%length(a)
% [vss,fval,exitFlag] = GlobalSolve(funcHd,n,-m*abs(max([L,h])),m*abs(max([L,h])));
% [vss,fval,exitFlag] = GlobalSolve(funcHd,n,0,m*abs(max([L,h])));
% a=1;
% ySolver(x)=vpa(yp(vss(1),vss(2),vss(3),vss(4),x))
% dySolver(x)=vpa(dyp(vss(1),vss(2),vss(3),vss(4),x));
% plotFunc(p1(1),p2(1),ySolver);
% hold on;
% plot([p1(1),p2(1)],[p1(2),p2(2)],'o');
% axis equal
% % plotFunc(p1(1),p2(1),dySolver);
% dy12=vpa([dyp(vss(1),vss(2),vss(3),vss(4),p1(1));dyp(vss(1),vss(2),vss(3),vss(4),p2(1))])
% exitFlag
% if output
%     chF=char(ySolver);
%     fileName=sprintf('GaussCurveL%dH%dA%dA%d.txt',L,h,angle1,angle2);
%     % fileName=sprintf('GaussCurveL%dH%dA%dA%.1f.txt',L,h,angle1,angle2);
%     fileID=fopen(fileName,'w');
%     fprintf(fileID,'%s\n',chF);
% end
% toc
% disp('!!!!!Done')