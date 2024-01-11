close all
clear 
clc
format long
size_m=400;
kappa=1;
%k=linspace(-pi,pi,size_m);
E=@(C,k) sqrt((2*C.*cos(k)).^2+kappa^2);

% index=0;
% C_x=0:0.01:3;
% for C=0:0.01:3
%     index=1+index;
%     E0(index)=E(C,0);
%     E1(index)=E(C,pi/4);
%     E2(index)=E(C,pi/3);
%     E3(index)=E(C,pi/2);
%     E4(index)=E(C,2*pi/3);
%     E5(index)=E(C,3*pi/4);
%     Epi(index)=E(C,pi);
% end
% figure(1);
% plot(C_x,E0)
% hold on
% plot(C_x,E1)
% hold on
% plot(C_x,E2)
% hold on
% plot(C_x,E3)
% hold on
% plot(C_x,E4)
% hold on
% plot(C_x,E5)
% hold on
% plot(C_x,Epi)
% 
% 
% plot(C_x,-E0)
% hold on
% plot(C_x,-E1)
% hold on
% plot(C_x,-E2)
% hold on
% plot(C_x,-E3)
% hold on
% plot(C_x,-E4)
% hold on
% plot(C_x,-E5)
% hold on
% plot(C_x,-Epi)
% xlabel("C");ylabel("E");legend("kappa=1")
% 

% index=0;
% C_x=0:0.001:0.1;
% for C=0:0.001:0.1
%     index=1+index;
%     E0(index)=E(C,0);
%     E1(index)=E(C,pi/4);
%     E2(index)=E(C,pi/3);
%     E3(index)=E(C,pi/2);
%     E4(index)=E(C,2*pi/3);
%     E5(index)=E(C,3*pi/4);
%     Epi(index)=E(C,pi);
% end
% figure(2);
% plot(C_x,E0)
% hold on
% plot(C_x,E1)
% hold on
% plot(C_x,E2)
% hold on
% plot(C_x,E3)
% hold on
% plot(C_x,E4)
% hold on
% plot(C_x,E5)
% hold on
% plot(C_x,Epi)
% 
% 
% plot(C_x,-E0)
% hold on
% plot(C_x,-E1)
% hold on
% plot(C_x,-E2)
% hold on
% plot(C_x,-E3)
% hold on
% plot(C_x,-E4)
% hold on
% plot(C_x,-E5)
% hold on
% plot(C_x,-Epi)
% xlabel("C");ylabel("E");legend("kappa=1")



index=0;
C_x=0:0.001:3;
for C=0:0.001:3
    index=1+index;
    E3(index)=E(C,pi/2);
end
figure(3);

plot(C_x,E3)

hold on
plot(C_x,-E3)
ylim([-3,3])
xlabel("C");ylabel("E");legend("kappa=1")