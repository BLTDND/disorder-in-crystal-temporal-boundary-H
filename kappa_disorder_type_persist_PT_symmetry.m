close all
clear
clc
format long
syms k
%% 保持手性
kappa=4;%2.3; 

k0=pi/3;%the central carrier wave frequency
C_before =4.4*1i;
C_after  = -7*1i;
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
%% 

%% constant

s=vg_before0./kappa/d_m/5;
%Lax-Friendrichs mathod
D0=(diag((1/2-s)*ones(N-1,1),1))+(diag((1/2+s)*ones(N,1),0));
% D0(end,1)=-s;

s=vg_after0*T/d_m/5;
D1=(diag((1/2-s)*ones(N/2-1,1),1))+(diag((1/2+s)*ones(N/2,1),0));
D1(end,1)=-s;

s=-s;
D2=(diag((1/2+s)*ones(N/2-1,1),-1))+(diag((1/2-s)*ones(N/2,1),0));
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
image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("constant")


%% periodical
figure(2)
for periodioc_level=1:12
    kappa=2+1.*cos(periodioc_level*k0.*(1:N));
    
    s=vg_before0./kappa/d_m/5;
    %Lax-Friendrichs mathod
    for i=1:N
        D0(i,i)=1/2+s(i);
        if i+1<N
            D0(i,i+1)=1/2-s(i);
        end
    end
    
    s=vg_after0./kappa/d_m/5;
    for i=1:N/2
        D1(i,i)=1/2+s(i);
        if i+1<N/2
            D1(i,i+1)=1/2-s(i);
        end
    end
    
    s=-s;
    for i=1:N/2
        D2(i,i)=1/2-s(i);
        if i-1>0
            D2(i,i-1)=1/2+s(i);
        end
    end
    
    
    
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
    subplot(3,4,periodioc_level);
    image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title(strcat("peridicial",num2str(periodioc_level)));
end


%% quasi-periodicity
kappa=2+fibonacci(N)./max(fibonacci(N));
s=vg_before0./kappa/d_m/5;
%Lax-Friendrichs mathod
for i=1:N
    D0(i,i)=1/2+s(i);
    if i+1<N
        D0(i,i+1)=1/2-s(i);
    end
end

s=vg_after0./kappa/d_m/5;
for i=1:N/2
    D1(i,i)=1/2+s(i);
    if i+1<N/2
        D1(i,i+1)=1/2-s(i);
    end
end

s=-s;
for i=1:N/2
    D2(i,i)=1/2-s(i);
    if i-1>0
        D2(i,i-1)=1/2+s(i);
    end
end



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
figure(3)
image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("qusi-peridicial")



%% aperiodicial
kappa=tmSequence(N, 2,6);
s=vg_before0./kappa/d_m/5;
%Lax-Friendrichs mathod
for i=1:N
    D0(i,i)=1/2+s(i);
    if i+1<N
        D0(i,i+1)=1/2-s(i);
    end
end

s=vg_after0./kappa/d_m/5;
for i=1:N/2
    D1(i,i)=1/2+s(i);
    if i+1<N/2
        D1(i,i+1)=1/2-s(i);
    end
end

s=-s;
for i=1:N/2
    D2(i,i)=1/2-s(i);
    if i-1>0
        D2(i,i-1)=1/2+s(i);
    end
end

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
figure(4)
image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("aperidocal")




%% disorder
figure(5)
psi_disorder=zeros(N/2,N,12);
for disorder_range=1:12
    for disorder_times=1:80
        kappa=4+(rand([1,N])-0.5)*8*(disorder_range-1)/12;
        s=vg_before0./kappa/d_m/5;
        for i=1:N
            D0(i,i)=1/2+s(i);
            if i+1<N
                D0(i,i+1)=1/2-s(i);
            end
        end
        
        s=vg_after0./kappa/d_m/5;
        for i=1:N/2
            D1(i,i)=1/2+s(i);
            if i+1<N/2
                D1(i,i+1)=1/2-s(i);
            end
        end
        
        s=-s;
        for i=1:N/2
            D2(i,i)=1/2-s(i);
            if i-1>0
                D2(i,i-1)=1/2+s(i);
            end
        end
        
        
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
        psi_disorder(:,:,disorder_range)=psi_disorder(:,:,disorder_range)+psi_2;
    end
    psi_disorder(:,:,disorder_range)=psi_disorder(:,:,disorder_range)./80;
    subplot(3,4,disorder_range),
    image(rescale(rot90(rot90(psi_disorder(1:2:end,:,disorder_range))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title(strcat("disorder=",num2str(disorder_range)))
end
% 
% image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("disorder")














