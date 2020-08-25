%%%%%%%%%%%%%%%%%%%%%%%%%
%y(x)=Gaussian Funtion
%y(0)=0,y'(0)=k1
%y(L)=-h0,y'(L)=k2
%%%%%%%%%%%%%%%%%%%%%%%%%
%	Grid caculate
%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% clc;
output=0;%output txt file or not
tic
syms y x h b c d 
y=h.*exp(-((x-b)./c).^2./2)+d; 
syms L h0 k1 k2
syms f1 f2
f1=h0+d.*(1-exp(-L.*(k2+k1.*(1+h0./d))./(2.*(h0+d))));
%%%Parameters set%%%%%%%%%%%%
h0_v=[31.62];
L_v=[156.85376343];
angle1=-5.19889199; k1_v=tand(angle1);
angle2=[-0.5]; k2_v=tand(angle2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lenH=length(h0_v);lenL=length(L_v);
lenAng1=length(angle1);lenAng2=length(angle2);
if ~(lenH==1 && lenL==1 && lenAng1==1 && lenAng2==1)
    [h0Grid,LGrid,k1Grid,k2Grid]=ndgrid(h0_v,L_v,k1_v,k2_v);
else
    fprintf('h0=%.4f\nL=%.4f\nk1=%.4f\nk2=%.4f\n',h0_v,L_v,k1_v,k2_v);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


f2=subs(f1,h0,h0Grid(:,1,1,1));
f2=subs(f2,L,LGrid(1,:,1,1));
f2=subs(f2,k1,k1Grid(1,1,:,1));
f1_d(d)=subs(f2,k2,k2Grid(1,1,1,:));
d_v=vpasolve(f1_d);
deltaD=vpa(f1_d(d_v));
c_v=sqrt(L_v.*(h0_v+d_v)./(k2_v-k1_v.*(1+h0_v./d_v)));
h_v=-d_v.*exp(c_v.^2.*k1_v.^2./(2.*d_v.^2));
b_v=-c_v.^2.*k1_v./d_v;
ySolver=subs(y,[h,b,c,d],[h_v,b_v,c_v,d_v]);
fprintf('\ndeltaD=%.4f\nh=%.4f\nb=%.4f\nc=%.4f\nd=%.4f\n',deltaD,h_v,b_v,c_v,d_v);
fprintf('\ny=%s\n',ySolver);
plotFunc(0,L_v,ySolver);
hold on;
plot([0,L_v],[0,-h0_v],'o');
axis equal

yS(x)=ySolver;y1(x)=diff(yS,x);
fprintf('\ny(0)=%.4f\ny(L)=%.4f\ny1(0)=%.4f\ny1(L)=%.4f\n',yS(0),yS(L_v),y1(0),y1(L_v));

if output
    chF=char(ySolver);
    fileName=sprintf('GaussCurveL%dH%dA%dA%d.txt',L_v,h0_v,angle1,angle2);
    fileID=fopen(fileName,'w');
    fprintf(fileID,'%s\n',chF);
end
toc
disp('!!!!!Done')