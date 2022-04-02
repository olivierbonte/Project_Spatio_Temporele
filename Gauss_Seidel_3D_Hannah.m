%% Gauss-Seidel 3D: modelleren aquifer
clear all
close all

%% Initialisaties
v_x = 0.1;
alpha_x = 0.5; 
D_x = alpha_x*v_x;
alpha_y = 0.05; 
D_y = alpha_y*v_x;
alpha_z = 0.05; %hetzelfde als alpha_y??
D_z = alpha_z*v_x;

% Onconditioneel stabiel volgens paper
Pe=0.5;
deltax=Pe*D_x/v_x;
deltay = deltax
deltaz = deltax;
Cr=0.5;
deltat = Cr*deltax/v_x;

x = 0:deltax:120; 
y = 0:deltay:14; 
z = 5:deltaz:15; %10m dikke aquifer startend op 5m diepte
t = 0:deltat:300; 
M = length(x)-2; 
N = length(y)-2;
O = length(z)-2;
K = length(t)-1;

u = zeros(N+4,M+4,O+4); %3D matrix, die telkens overschreden zal worden

% Initiële conditie
idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
idz = z <= 7; %in de bovenste 2m
posx = find(idx) + 1; %rekening houdend met imaginaire knopen
posy = find(idy) + 1;
posz = find(idz) + 1;
% Balk met C0=1 als initiële conditie?
C0 = 1; %bv 1 kg/m^3 dus
u(posy,posx,posz) = C0;

%Of bijv. lineaire afname tussen 5 en 7m diep 
for g=posz
    u(posy,posx,g)=1-0.5*(z(g-1)-5);
end

%% Plot initiele conditie
% figure()
% surf(x,y,z,u(2:end-1,2:end-1,2:end-1,1),'LineStyle','none')
% figure()
% contourf(x,y,u(2:end-1,2:end-1,1))
% xlabel('x')
% ylabel('y')
% colorbar

%% Gauss Seidel iteratie
difference = 1;
phi = 1e-4; %convergentiecriterium
iter = 0;

% Coeffiënten
R = 1; % Geen retardatie
%LL
B1 = R/deltat + 2*D_x/deltax^2 + 2*D_y/deltay^2 + 2*D_z/deltaz^2;
%RL
B2 = D_x/deltax^2 - v_x/(2*deltax);
B3 = D_x/deltax^2 + v_x/(2*deltax);
B4 = D_y/deltay^2;
B5 = B4;
B6 = D_z/deltaz^2;
B7 = B6;

unew = u;
% u_flat = u(:,:,2); %laat deze dus op de initiële condities starten met z'n iteraties!
% unew_flat = u_flat;
% figure()
for k = 1:K %tot tijdstip K+1 door te indexeren op k+1
    iter = 0;
    while difference > phi
        for g = 2:O+3 %van g=0 tot g=O+1 
            for j = 2:N+3 %van j=0 tot j=N+1
                for i = 2:M+3 %van i=0 tot i=M+1
                    unew(j,i,g) = B2/B1*u(j,i+1,g) + B3/B1*unew(j,i-1,g) +...
                        B4/B1*u(j+1,i,g) + B5/B1*unew(j-1,i,g) + B6/B1*u(j,i,g+1) +...
                        B7/B1*unew(j,i,g-1) + R/(deltat*B1)*u(j,i,k);
                end
                unew_flat(j,1) = unew_flat(j,3); %Neumann links: u_-1 = u_1
                unew_flat(j,M+4) = unew_flat(j,M+2); %Neuman rehcts: u_M+2 = u_M
            end
            unew_flat(1,:) = unew_flat(3,:); %Neumann onderaan
            unew_flat(N+4,:) = unew_flat(N+2,:); %Neumann bovenaan
            difference = max(max(abs(unew_flat-u_flat)));
            u_flat = unew_flat;
            iter = iter +1;
        end
    end
    iter;
    k;
    u(:,:,k+1) = u_flat;
    %contourf(x,y,u(2:end-1,2:end-1,k+1));
    %xlabel('x')
    %ylabel('y')
    %zlabel('u')
    %title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    %colorbar
    %Ani(k) = getframe;
    if rem(k,50) == 0
       k 
    end
    difference = 1; %zodat hij terug vertrekt
end
%% Opslaan

%% visualisaties
%videostijl
figure()
for k = 0:10:K %Dus niet alle tijdstappen visualiseren
    contourf(x,y,u(2:end-1,2:end-1,k+1));
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    colorbar
    Ani(k+1) = getframe;
end

tijden = [0,200,600];
f = figure();
f.Position(3:4) = [1.5*560,1.2*420];
for i = 1:length(tijden)
    tijdid = t == tijden(i);
    subplot(3,1,i)
    contourf(x,y,u(2:end-1,2:end-1,tijdid))
    xlabel('x')
    ylabel('y')
    title(strcat('t = ',num2str(t(tijdid))))
    colorbar
end
exportgraphics(gcf,'Figuur_9.png','Resolution',900)