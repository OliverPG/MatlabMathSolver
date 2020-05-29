clear;clc;
output=0;%output txt file
tic
L=175./1000;
h=150./1000;
p1=[0,0]; angle1=0;
k1=tand(angle1);
p2=[L,-h]; angle2=2;
k2=-tand(angle2);
syms y x a b c d
% y(x)=a*exp(-(x-b)^2/(2*c^2))+d;
y(x)=a*exp(-c*(x-b)^2)+d;
% y(x)=a*exp(-(x*cos(x)-b)^2/(2*c^2))+d; 
yp(a,b,c,d,x)=y(x);
dy(x)=diff(y,x); dyp(a,b,c,d,x)=dy(x);
func1(a,b,c,d)=y(p1(1))-p1(2);
func2(a,b,c,d)=y(p2(1))-p2(2);
func3(a,b,c,d)=dy(p1(1))-k1;
func4(a,b,c,d)=dy(p2(1))-k2;
funcs=[func1(a,b,c,d);func2(a,b,c,d);func3(a,b,c,d);func4(a,b,c,d)];
vs=[a,b,c,d];
vs0=[h,L,1.5*L/2,h];
funcHd=matlabFunction(funcs,'vars',{[a,b,c,d]});
n=4;m=20;
[vss,fval,exitFlag] = GlobalSolve(funcHd,n);
% [vss,fval,exitFlag] = GlobalSolve(funcHd,n,-m*abs(max([L,h])),m*abs(max([L,h])));
% [vss,fval,exitFlag] = GlobalSolve(funcHd,n,0,m*abs(max([L,h])));
ySolver(x)=vpa(yp(vss(1),vss(2),vss(3),vss(4),x))
dySolver(x)=vpa(dyp(vss(1),vss(2),vss(3),vss(4),x));
plotFunc(p1(1),p2(1),ySolver);
hold on;
plot([p1(1),p2(1)],[p1(2),p2(2)],'o');
axis equal
% plotFunc(p1(1),p2(1),dySolver);
dy12=vpa([dyp(vss(1),vss(2),vss(3),vss(4),p1(1));dyp(vss(1),vss(2),vss(3),vss(4),p2(1))])
exitFlag
if output
    chF=char(ySolver);
    fileName=sprintf('GaussCurveL%dH%dA%dA%d.txt',L,h,angle1,angle2);
    % fileName=sprintf('GaussCurveL%dH%dA%dA%.1f.txt',L,h,angle1,angle2);
    fileID=fopen(fileName,'w');
    fprintf(fileID,'%s\n',chF);
end
toc
disp('!!!!!Done')