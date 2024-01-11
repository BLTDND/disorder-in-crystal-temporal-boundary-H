close all
clear
clc
format long
syms k
kappa=2.3; 

k0=pi/3;%the central carrier wave frequency
C_before =4.4;
C_after  = -7;
C=C_before;
E=sqrt((2*C*cos(k))^2+kappa^2);
vg_before=diff(E,k);
vg_before0=double(subs(vg_before,k,pi/3));
C=C_after;
E=sqrt((2*C*cos(k))^2+kappa^2);
vg_after=diff(E,k);
vg_after0=double(subs(vg_after,k,pi/3));% the group velocity

N=200;
m=linspace(-200,200,N);%for standing wave
d_m=m(2)-m(1);
m_value0=N/2*d_m;
sigma=23;

T=1/kappa;

%% initial state 
phi_plus=@(k,C) [2*C.*cos(k)+sqrt((2*C.*cos(k)).^2+kappa^2); kappa]./sqrt((2*C.*cos(k)+sqrt((2*C.*cos(k)).^2+kappa^2)).^2+kappa^2);
phi_subtract=@(k,C) [2*C.*cos(k)-sqrt((2*C.*cos(k)).^2+kappa^2); kappa]./sqrt((2*C.*cos(k)-sqrt((2*C.*cos(k)).^2+kappa^2)).^2+kappa^2);
phi_initial=phi_plus(k0,C_before);
phi_0=double(phi_initial);
phi_0=[1;1]/sqrt(2);
N_t=100;

E_initial=zeros(1,N);
for m_index=1:N
    m_value=m_index*d_m;
    if mod(m_index,2)==1
        E_initial(1,m_index)=1/sqrt(sigma*sqrt(pi))*exp(-(m_value-m_value0)^2/2/sigma^2)*exp(1i*k0*m_value).*double(phi_initial(1)); % Eqn.21 
    else
        E_initial(1,m_index)=1/sqrt(sigma*sqrt(pi))*exp(-(m_value-m_value0)^2/2/sigma^2)*exp(1i*k0*m_value).*double(phi_initial(2)); % Eqn.21 
    end
end
E_initial=norm_matrix(E_initial);
%% ini
s=vg_before0*T/d_m;
%Lax-Friendrichs mathod
D0=(diag((-s)*ones(N-1,1),1))+(diag((1+s)*ones(N,1),0));
% D0(end,1)=-s;

s=vg_after0*T/d_m;
D1=(diag((-s)*ones(N/2-1,1),1))+(diag((1+s)*ones(N/2,1),0));
D1(end,1)=-s;

s=-s;
D2=(diag((s)*ones(N/2-1,1),-1))+(diag((1-s)*ones(N/2,1),0));
D2(1,end)=s;


psi=zeros(N_t,N);
psi(1,:)=E_initial;
for T_index=2:N_t 
    t=T*(T_index-N_t/2);
    if T_index<N_t/2
        psi(T_index,:)=psi(T_index-1,:)*D0;
    else
        psi_odd=psi(T_index-1,1:2:end);
        psi_even=psi(T_index-1,2:2:end);
        psi(T_index,1:2:end)=psi_odd*D1;
        psi(T_index,2:2:end)=psi_even*D2;
%         psi(:,T_index)=D2*psi(:,T_index-1)+D1*psi(:,T_index-1);
    end
    psi(T_index,:)=norm_matrix(psi(T_index,:));
end
psi_2=abs(psi).^2;
figure(1)
image(rescale(rot90(rot90(psi_2)),0,255));xlabel('omiga_m');ylabel('t');colorbar


























