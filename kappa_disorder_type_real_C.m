close all
clear
clc
format long
syms k
% real all 准确相移
kappa=6;%2.3; 

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
% phi_0=[1;1]/sqrt(2);
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

s=vg_before0*T/d_m/5;
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
for periodic_level=1:12
    kappa=20+6.*cos(periodic_level*k0.*(1:N));
    
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
    psi_periodic(:,:,periodic_level)=psi_2;
    subplot(4,3,periodic_level);
    image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title(strcat("periodicial",num2str(periodic_level)));
end

for disorder_range=1:12
    tr_before=find(psi_periodic(50,:,disorder_range)==max(psi_periodic(50,:,disorder_range)));
    Incident_angle(disorder_range)=pi/2-atan(abs(tr_before-200)/50);
    tr_after=find(psi_periodic(end,:,disorder_range)==max(psi_periodic(end,tr_before:end,disorder_range)));
    Reflection_angle(disorder_range)=pi/2-atan(abs(tr_after-tr_before)/50);
    re_after=find(psi_periodic(end,:,disorder_range)==max(psi_periodic(end,1:tr_before,disorder_range)));
    Refraction_angle(disorder_range)=atan(abs(re_after-tr_before)/50);
end
level=1:12;
order =4;
coefficients_Incident_angle = polyfit(level, Incident_angle, order);

% 绘制拟合曲线1
xFit_Incident_angle= linspace(min(level), max(level), 100);
yFit_Incident_angle = polyval(coefficients_Incident_angle, xFit_Incident_angle);


order =4;
coefficients_Reflection_angle= polyfit(level, Reflection_angle, order);

% 绘制拟合曲线2
xFit_Reflection_angle= linspace(min(level), max(level), 100);
yFit_Reflection_angle= polyval(coefficients_Reflection_angle, xFit_Reflection_angle);



order =4;
coefficients_Refraction_angle= polyfit(level, Refraction_angle, order);

% 绘制拟合曲线3
xFit_Refraction_angle= linspace(min(level), max(level), 100);
yFit_Refraction_angle = polyval(coefficients_Refraction_angle, xFit_Refraction_angle);



figure(11)


plot(Incident_angle,"*r")

hold on
plot(Reflection_angle,"ob")

hold on
hold on
plot(Refraction_angle,"+g")





plot(xFit_Incident_angle,yFit_Incident_angle,"-r")
hold on
plot(xFit_Reflection_angle,yFit_Reflection_angle,"-b")
hold on
plot(xFit_Refraction_angle,yFit_Refraction_angle,"-g")% 
% image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("disorder")
legend("Incident angle","Reflection angle","Refraction angle","FIT Incident angle","FIT Reflection angle","FIT Refraction angle")
xlabel("level")
ylabel("angle")

%% quasi-periodicity
kappa=14+fibonacci(N);
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
image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("qusi-periodicial")



%% aperiodicial
kappa=tmSequence(N, 14,20);
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
image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("aperiodical")




%% disorder
figure(5)
psi_disorder=zeros(N/2,N,12);
for disorder_range=1:12
    for disorder_times=1:80
        kappa=20+(rand([1,N])-0.5)*10*(disorder_range-1)/12;
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


for disorder_range=1:12
    tr_before=find(psi_disorder(50,:,disorder_range)==max(psi_disorder(50,:,disorder_range)));
    Incident_angle(disorder_range)=pi/2-atan(abs(tr_before-200)/50);
    tr_after=find(psi_disorder(end,:,disorder_range)==max(psi_disorder(end,tr_before:end,disorder_range)));
    Reflection_angle(disorder_range)=pi/2-atan(abs(tr_after-tr_before)/50);
    re_after=find(psi_disorder(end,:,disorder_range)==max(psi_disorder(end,1:tr_before,disorder_range)));
    Refraction_angle(disorder_range)=atan(abs(re_after-tr_before)/50);
end
level=1:12;
order =4;
coefficients_Incident_angle = polyfit(level, Incident_angle, order);

% 绘制拟合曲线1
xFit_Incident_angle= linspace(min(level), max(level), 100);
yFit_Incident_angle = polyval(coefficients_Incident_angle, xFit_Incident_angle);


order =4;
coefficients_Reflection_angle= polyfit(level, Reflection_angle, order);

% 绘制拟合曲线2
xFit_Reflection_angle= linspace(min(level), max(level), 100);
yFit_Reflection_angle= polyval(coefficients_Reflection_angle, xFit_Reflection_angle);



order =4;
coefficients_Refraction_angle= polyfit(level, Refraction_angle, order);

% 绘制拟合曲线3
xFit_Refraction_angle= linspace(min(level), max(level), 100);
yFit_Refraction_angle = polyval(coefficients_Refraction_angle, xFit_Refraction_angle);



figure(6)


plot(Incident_angle,"*r")

hold on
plot(Reflection_angle,"ob")

hold on
hold on
plot(Refraction_angle,"+g")





plot(xFit_Incident_angle,yFit_Incident_angle,"-r")
hold on
plot(xFit_Reflection_angle,yFit_Reflection_angle,"-b")
hold on
plot(xFit_Refraction_angle,yFit_Refraction_angle,"-g")% 
% image(rescale(rot90(rot90(psi_2(1:2:end,:))),0,255));xlabel('omiga_m');ylabel('t');colorbar,title("disorder")
legend("Incident angle","Reflection angle","Refraction angle","FIT Incident angle","FIT Reflection angle","FIT Refraction angle")
xlabel("level")
ylabel("angle")













