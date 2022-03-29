%% poging to expliciet op te lossen!
clear all
close all
tic

%Dit faalt: volgens stabilteitsvereiste van de cursus p 62 wordt de deltat
%veel te klein => gigantisch aantal tijdstappen nodig... 
%krijg wel een plot voor kortere delta t 
%% Initialisaties
v = 0.1; %m/d
alpha_L = 0.05; %m
D_L = alpha_L*v;
Pe = 1; %Peclet nummer
deltax = Pe*alpha_L;

%Bepaal hier deltat ahv delta = 1/2*deltax^2
deltat = (1/2-1e-5)*deltax^2;

alpha_T = 0.005; %m
D_T = alpha_T*v;

%grid
x = 0:deltax:120; %x tot 120 m
%bij gebrek aan info: neem hier deltay = deltax
deltay = deltax;
y = 0:deltay:14; %y tot 14m
t = 0:deltat:10; %tot dag 600 in de figuur: GAAT niet
M = length(x)-2; %aftrekken van interne knopen
N = length(y)-2;
K = length(t)-1;

R =1; %dus geen extra retardatie

u = zeros(N+2,M+2,K+1);

idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
C0 = 1; %bv 1 kg/m^3 dus
u(idy,idx,1) = C0; %de initiële condities

figure()
for k = 1:K
    for j = 1:N+2
        for i = 1:M+2
            if i == 1 & j == 1
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i+1,k) -2*u(j,i,k) + u(j,i+1,k)) + ...
                    D_T/deltay^2*(u(j+1,i,k) - 2*u(j,i,k) + u(j+1,i,k)) - ...
                    v/(2*deltax)*(u(j,i+1,k)-u(j,i+1,k))) + u(j,i,k);
            elseif i ==1 & j == N+2
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i+1,k) -2*u(j,i,k) + u(j,i+1,k)) + ...
                    D_T/deltay^2*(u(j-1,i,k) - 2*u(j,i,k) + u(j-1,i,k)) - ...
                    v/(2*deltax)*(u(j,i+1,k)-u(j,i+1,k))) + u(j,i,k);
            elseif i == 1 
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i+1,k) -2*u(j,i,k) + u(j,i+1,k)) + ...
                    D_T/deltay^2*(u(j+1,i,k) - 2*u(j,i,k) + u(j-1,i,k)) - ...
                    v/(2*deltax)*(u(j,i+1,k)-u(j,i+1,k))) + u(j,i,k);
            elseif i == M+2 & j == 1
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i-1,k) -2*u(j,i,k) + u(j,i-1,k)) + ...
                    D_T/deltay^2*(u(j+1,i,k) - 2*u(j,i,k) + u(j+1,i,k)) - ...
                    v/(2*deltax)*(u(j,i-1,k)-u(j,i-1,k))) + u(j,i,k);
            elseif i == M+2 & j == N+2
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i-1,k) -2*u(j,i,k) + u(j,i-1,k)) + ...
                    D_T/deltay^2*(u(j-1,i,k) - 2*u(j,i,k) + u(j-1,i,k)) - ...
                    v/(2*deltax)*(u(j,i-1,k)-u(j,i-1,k))) + u(j,i,k);
            elseif i == M+2
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i-1,k) -2*u(j,i,k) + u(j,i-1,k)) + ...
                    D_T/deltay^2*(u(j+1,i,k) - 2*u(j,i,k) + u(j-1,i,k)) - ...
                    v/(2*deltax)*(u(j,i-1,k)-u(j,i-1,k))) + u(j,i,k);
            elseif j == 1 
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i+1,k) -2*u(j,i,k) + u(j,i-1,k)) + ...
                    D_T/deltay^2*(u(j+1,i,k) - 2*u(j,i,k) + u(j+1,i,k)) - ...
                    v/(2*deltax)*(u(j,i+1,k)-u(j,i-1,k))) + u(j,i,k);
            elseif j == N+2
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i+1,k) -2*u(j,i,k) + u(j,i-1,k)) + ...
                    D_T/deltay^2*(u(j-1,i,k) - 2*u(j,i,k) + u(j-1,i,k)) - ...
                    v/(2*deltax)*(u(j,i+1,k)-u(j,i-1,k)));
            else
                u(j,i,k+1) = deltat/R*(D_L/deltax^2*(u(j,i+1,k) -2*u(j,i,k) + u(j,i-1,k)) + ...
                    D_T/deltay^2*(u(j+1,i,k) - 2*u(j,i,k) + u(j-1,i,k)) - ...
                    v/(2*deltax)*(u(j,i+1,k)-u(j,i-1,k))) + u(j,i,k);
            end
        end
    end
    contourf(x,y,u(:,:,k+1));
    colorbar
    xlabel('x')
    ylabel('y')
    zlabel('u')
    colorbar
    title(strcat('t = ',num2str(t(k))))
    shg
    pause(deltat)
end



tijden = 1:10;
for i = 1:length(tijden)
    pos = find(t >= tijden(i));
    tijdid = pos(1);
    figure()
    contourf(x,y,u(:,:,tijdid))
    colorbar
    xlabel('x')
    ylabel('y')
    title(strcat('t =',num2str(t(tijdid))))
end
toc