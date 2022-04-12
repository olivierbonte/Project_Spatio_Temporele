%% Figuur 9 via Gauss-Seidel
%clear all
close all
tic 
%aanpassing naar 2D matrix omdat anders te groot! 
%% Initialisaties
v = 0.1; %m/d
alpha_L = 0.05; %m
D_L = alpha_L*v;
Pe = 0.5; %Peclet nummer nemen op 0.5 zodat kleiner dan 1. 
deltax = Pe*alpha_L; %deltax = Pe*D_L/v = Pe*alpha_L

Cr = 0.5; %Courant nummer zodat kleiner dan 1!
deltat = Cr*deltax/v; %in dagen

alpha_T = 0.005; %m
D_T = alpha_T*v;

%grid
x = 0:deltax:120; %x tot 120 m
%bij gebrek aan info: neem hier deltay = deltax
deltay = deltax;
y = 0:deltay:14; %y tot 14m
t = 0:deltat:600; %tot dag 600 in de figuur
M = length(x)-2; %aftrekken van interne knopen
N = length(y)-2;
K = length(t)-1;

u = zeros(N+4,M+4); %geen 3D matrix
%initiele condities
idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
%beter idee! om rekening te houden met de imaginaire knopen
posx = find(idx) + 1; %omdat indexering start op -1
posy = find(idy) +1; %idem
C0 = 1; %bv 1 kg/m^3 dus
u(posy,posx) = C0; %Dit is tijdstip 1
uprevious = u; %Dit is tijdstip 0
u0 = u; %opslaan voor later! 

%plot de initiele conditie
figure()
surf(x,y,u(2:end-1,2:end-1,1),'LineStyle','none')
figure()
contourf(x,y,u(2:end-1,2:end-1,1))
xlabel('x')
ylabel('y')
colorbar

%% Gauss Seidel iteratie
difference = 1;
phi = 1e-4; %convergentiecriterium
iter = 0;
%preallocaties
R = 1; %dus geen retardatie hier
B1 = R/deltat + 2*D_L/deltax^2 + 2*D_T/deltay^2;
B2 = -D_L/deltax^2 + v/(2*deltax);
B3 = -D_L/deltax^2 - v/(2*deltax);
B4 = -D_T/deltay^2;
B5 = B4;

unew = u; %belangrijke preallocatie
figure()
uani = {};
ani_index = 1;
for k = 1:K %tot tijdstip K (K+1) maar indexeren op k+1
    iter = 0;
    while difference > phi
        for j = 2:N+3 % van j = 0 tot j = N+1
            for i = 2:M+3 % van i = 0 tot i = M+1
                unew(j,i) = -B2/B1*unew(j,i+1)-B3/B1*unew(j,i-1)-...
                    B4/B1*unew(j+1,i) - B5/B1*unew(j-1,i) + R/(deltat*B1)*uprevious(j,i);
            end
            unew(j,1) = unew(j,3); %Neumann links: u_-1 = u_1
            unew(j,M+4) = unew(j,M+2); %Neuman rehcts: u_M+2 = u_M
        end
        unew(1,:) = unew(3,:); %Neumann onderaan
        unew(N+4,:) = unew(N+2,:); %Neumann bovenaan
        difference = max(max(abs(unew-u)));
        u = unew;
        iter = iter +1;
    end        
    uprevious = u; %belang om hier op te slaan!
    if rem(t(k+1),10) == 0
        k
        contourf(x,y,u(2:end-1,2:end-1));
        xlabel('x')
        ylabel('y')
        zlabel('u')
        title(strcat('tijd =  ',num2str(t(k+1)),'dagen'))
        colorbar
        Ani(k+1) = getframe;
        uani{ani_index} = u;
        ani_index = ani_index + 1;
    end
    if t(k+1)>=200-deltat/2 && t(k+1)<=200+deltat/2
        u200 = u;
    elseif t(k+1)>=600-deltat
        u600 = u;
    end
    difference = 1; %zodat hij terug vertrekt
end
%% Opslaan

%% visualisaties
%videostijl
figure()
lengte = length(uani);
for m = 1:lengte %Dus niet alle tijdstappen visualiseren
    u = uani{m};
    contourf(x,y,u(2:end-1,2:end-1));
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    colorbar
    Ani(k+1) = getframe;
end
%% plot voor rapport 
tijden = [0,200,600];
f = figure();
f.Position(3:4) = [1.5*560,1.2*420];
us = {u0,u200,u600};
for i = 1:length(us)
    %tijdid = t == tijden(i);
    subplot(3,1,i)
    u = us{i};
    contourf(x,y,u(2:end-1,2:end-1))
    xlabel('x [m]')
    ylabel('y [m]')
    title(strcat('t = ', num2str(tijden(i)) ,' days'))
    colorbar
end
exportgraphics(gcf,'Figuur_9.png','Resolution',900)
toc