close all 
clear all

%% Stabiel

deltax=5e-2;
deltat=1.2e-3;
deltat/deltax^2 %dus KLEINER dan 0.5 en stabiel
%Omega=; %niet nodig hier?
Time=41*deltat; %tot 40s dus
x=0:deltax:1;
t=0:deltat:Time;
M=length(x)-2; %is dus aantal INTERNE knopen
K=length(t)-1; %is dus aantal knopen -1 want indexering start op 0
u=zeros(K+1,M+2);


u(1,:) = 4.*x.*(1-x); %initiële
u(2:end,M+2) = 0; %Dirichlet. geldt pas vanaf t > 0
u(2:end,1) = 0;



for j= 1:K %Want initieel gekend. Start op j omdat werken met j+1 index!
    %gaan tem index K => K+1 in matlab
    for i=2:M+1 % van i = 1 tem i = M. 0 en M+1 gekend
        u(j+1,i) = u(j,i) + deltat/deltax^2*(u(j,i+1)-2*u(j,i)+u(j,i-1));
    end
end

figure()
surf(x,t,u,'linestyle','none')
%surf(x,t,u(:,2:M+3),'linestyle','none')
title(strcat('\Delta t =', num2str(deltat)))
xlabel('x')
ylabel('t')
zlabel('u')

figure()
indices = [1,21,41];
plot(x,u(indices,:))
legend('t^0','t^{20}','t^{40}')

%% ONstabiel

deltax=5e-2;
deltat=1.4e-3;
deltat/deltax^2 %dus groter dan 0.5 en NIET stabiel
%Omega=; %niet nodig hier?
Time=41*deltat; %tot 40s dus
x=0:deltax:1;
t=0:deltat:Time;
M=length(x)-2; %is dus aantal INTERNE knopen
K=length(t)-1; %is dus aantal knopen -1 want indexering start op 0
u=zeros(K+1,M+2);


u(1,:) = 4.*x.*(1-x); %initiële
u(2:end,M+2) = 0; %Dirichlet. geldt pas vanaf t > 0
u(2:end,1) = 0;



for j= 1:K %Want initieel gekend. Start op j omdat werken met j+1 index!
    %gaan tem index K => K+1 in matlab
    for i=2:M+1 % van i = 1 tem i = M. 0 en M+1 gekend
        u(j+1,i) = u(j,i) + deltat/deltax^2*(u(j,i+1)-2*u(j,i)+u(j,i-1));
    end
end

figure()
surf(x,t,u,'linestyle','none')
%surf(x,t,u(:,2:M+3),'linestyle','none')
title(strcat('\Delta t =', num2str(deltat)))
xlabel('x')
ylabel('t')
zlabel('u')

figure()
indices = [1,21,41];
plot(x,u(indices,:))
legend('t^0','t^{20}','t^{40}')