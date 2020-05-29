clear;clc;
L=1473;
h=91.27;
p1=[0,0]; k1=tand(30);
p2=[L,-h]; k2=-tand(30);
syms y x a b c d
y(x)=a*exp(-(x-b)^2/(2*c^2))+d; yp(a,b,c,d,x)=y(x);
dy(x)=diff(y,x); dyp(a,b,c,d,x)=dy(x);
func1(a,b,c,d)=y(p1(1))-p1(2);
func2(a,b,c,d)=y(p2(1))-p2(2);
func3(a,b,c,d)=dy(p1(1))-k1;
func4(a,b,c,d)=dy(p2(1))-k2;
funcs=[func1(a,b,c,d);func2(a,b,c,d);func3(a,b,c,d);func4(a,b,c,d)];
vs=[a,b,c,d];
vs0=[h,L,1.5*L/2,h];
% vs0=[h,0,1,-h];
funcHd=matlabFunction(funcs,'vars',{[a,b,c,d]});
tolerance=1e-16;
tolFunc=1e-16;
tolStep=1e-16;
maxCacuN=3e4;
options = optimoptions('fsolve','Display','none','Algorithm','levenberg-marquardt','CheckGradients',true,'FunctionTolerance',tolFunc,'OptimalityTolerance',tolerance,'StepTolerance',tolStep,'Maxfuneval',maxCacuN,'MaxIterations',maxCacuN,'PlotFcn',@optimplotfirstorderopt);
% options = optimoptions('fsolve','Display','none','Algorithm','trust-region-dogleg','FunctionTolerance',tolerance,'OptimalityTolerance',tolerance,'StepTolerance',tolerance,'Maxfuneval',20000,'MaxIterations',20000,'PlotFcn',@optimplotfirstorderopt);
% options = optimoptions('fsolve','Display','none','Algorithm','trust-region','FunctionTolerance',tolerance,'OptimalityTolerance',tolerance,'StepTolerance',tolerance,'Maxfuneval',20000,'MaxIterations',20000,'PlotFcn',@optimplotfirstorderopt);
% options = optimoptions('fsolve','Display','none','Algorithm','trust-region-reflective','FunctionTolerance',tolerance,'OptimalityTolerance',tolerance,'StepTolerance',tolerance,'Maxfuneval',20000,'MaxIterations',20000,'PlotFcn',@optimplotfirstorderopt);
% 'levenberg-marquardt' 'trust-region-dogleg' 'trust-region'
options = optimset('GradObj','on');
options = optimset('Hessian','on')
% options.MaxFunctionEvaluations=4000;
% [vss,fval,exitflag,output]=fsolve(funcHd,vs0,options);
[vss,fval,exitflag,output]=fsolve(funcHd,vs0);
% [vss,fval,exitflag,output]=fminsearch(funcHd,vs0,options)
ySolver(x)=vpa(yp(vss(1),vss(2),vss(3),vss(4),x))
plotFunc(p1(1),p2(1),ySolver);
hold on;
plot([p1(1),p2(1)],[p1(2),p2(2)],'o');
axis equal
dy12=vpa([dyp(vss(1),vss(2),vss(3),vss(4),p1(1));dyp(vss(1),vss(2),vss(3),vss(4),p2(1))])
char(ySolver)