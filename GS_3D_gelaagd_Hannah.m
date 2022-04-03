%% Gauss-Seidel 3D: modelleren gelaagde aquifer
clear all
close all

%% Initialisaties
% Bovenste laag: loamy sand (2-11m)
% Onderste laag: silt loam (11-20m)

K_hydr=[12.4*10^(-6), 3.9*10^(-6)]; %hydraulische geleidbaarheid in m/s
K_hydr=K_hydr*60*60*24; %m/d
v_x_2=[0.1 0.1*K_hydr(2)/K_hydr(1)];
% Veronderstel de oorspronkelijke snelheid in de bovenste laag.
% De snelheid in de onderste laag wordt dan berekend o.b.v. verhoudingen.

alpha=[0.001 0.067];
alpha_x=[0.2 0.2*alpha(2)/alpha(1)]; %m
% Opm: 500 keer groter dan in de paper, vanwege beperkte opslag voor 3D
% matrix...
D_x_2=alpha_x.*v_x_2;
alpha_y=alpha_x/10; %alpha_y en alpha_z opnieuw 10x kleiner??
D_y_2=alpha_y.*v_x_2;
D_z_2=D_y_2;

Pe=0.8; %zie opm
deltax=round(Pe*min(alpha_x),2); %kleinere van de 2 alpha_x gebruiken
deltay=deltax;
deltaz=deltax;
Cr=0.8; %zie opm
deltat=round(Cr*deltax/max(v_x_2),2); %grotere van de 2 v_x gebruiken

x = 0:deltax:120; 
y = 0:deltay:14; 
z = 2:deltaz:20; %18m dikke aquifer startend op 2m diepte (GWT)
t = 0:deltat:600; 
M = length(x)-2; 
N = length(y)-2;
O = length(z)-2;
K = length(t)-1;

u = zeros(N+4,M+4,O+4); %3D matrix, die telkens overschreven zal worden

% Initiële conditie: balk
idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
idz = z >= 9 & z <= 13; 
posx = find(idx) + 1; %rekening houdend met imaginaire knopen
posy = find(idy) + 1;
posz = find(idz) + 1;
% Balk met C0=1 als initiële conditie?
C0 = 1; %bv 1 kg/m^3 dus
u(posy,posx,posz) = C0;

%% Plot initiele conditie
% In xz vlak
id=find(y>(7-deltay/2) & y<(7+deltay/2))+1;
u_xz = squeeze(u(id(1),2:end-1,2:end-1)); %om om te zetten naar 2D matrix
figure(1)
contourf(x,z,u_xz') %transponeren van u_xz!
xlabel('x')
ylabel('z')
colorbar
set (gca, 'ydir', 'reverse') %omdraaien van y-as!
title('Initial conditions in the xz-plane: y = 7 m')

%% Gauss Seidel iteratie
difference = 1;
phi = 1e-4; %convergentiecriterium
iter = 0;

R = 1; % Geen retardatie

u_k = u; %3D matrix van de vorige tijdstap (k) initialiseren
unew = u; %3D matrices van de huidige tijdstap (k+1) initialiseren

figure(2)
hold on
for k = 1:K %tot tijdstip K+1 door te indexeren op k+1
    iter = 0;
    while difference > phi
        for g = 2:O+3 %van g=0 tot g=O+1 
            if z(g-1)<=11
                % v en D
                v_x=v_x_2(1);
                D_x=D_x_2(1);
                D_y=D_y_2(1);
                D_z=D_z_2(1);
                % Coefficienten
                B1 = R/deltat + 2*D_x/deltax^2 + 2*D_y/deltay^2 + 2*D_z/deltaz^2;
                B2 = D_x/deltax^2 - v_x/(2*deltax);
                B3 = D_x/deltax^2 + v_x/(2*deltax);
                B4 = D_y/deltay^2;
                B5 = B4;
                B6 = D_z/deltaz^2;
                B7 = B6;
            elseif z(g-1)>11
                % v en D
                v_x=v_x_2(2);
                D_x=D_x_2(2);
                D_y=D_y_2(2);
                D_z=D_z_2(2);
                % Coefficienten
                B1 = R/deltat + 2*D_x/deltax^2 + 2*D_y/deltay^2 + 2*D_z/deltaz^2;
                B2 = D_x/deltax^2 - v_x/(2*deltax);
                B3 = D_x/deltax^2 + v_x/(2*deltax);
                B4 = D_y/deltay^2;
                B5 = B4;
                B6 = D_z/deltaz^2;
                B7 = B6;
            end
            for j = 2:N+3 %van j=0 tot j=N+1
                for i = 2:M+3 %van i=0 tot i=M+1
                    unew(j,i,g) = B2/B1*u(j,i+1,g) + B3/B1*unew(j,i-1,g) +...
                        B4/B1*u(j+1,i,g) + B5/B1*unew(j-1,i,g) + B6/B1*u(j,i,g+1) +...
                        B7/B1*unew(j,i,g-1) + R/(deltat*B1)*u_k(j,i,g);
                end
                unew(j,1,g) = unew(j,3,g); %von Neumann randvoorwaarden naar x
                unew(j,M+4,g) = unew(j,M+2,g); 
            end
            unew(1,:,g) = unew(3,:,g); %von Neumann randvoorwaarden naar y
            unew(N+4,:,g) = unew(N+2,:,g); 
        end
        unew(:,:,1)=unew(:,:,3); %von Neumann randvoorwaarden naar z
        unew(:,:,O+4)=unew(:,:,O+2);
        difference = max(max(max(abs(unew-u))));
        u = unew;
        iter = iter +1;
    end
    iter;
    k;
    u_k=u;

    % Animatie figuur
    if rem(k,20)==0
        id=find(y>(7-deltay/2) & y<(7+deltay/2))+1;
        u_xz = squeeze(u(id(1),2:end-1,2:end-1)); %om om te zetten naar 2D matrix
        contourf(x,z,u_xz') %transponeren van u_xz!
        xlabel('x')
        ylabel('z')
        set (gca, 'ydir', 'reverse')
        colorbar
        title(strcat('GS in xz-plane: y = 7 m, t = ',num2str(t(k+1)),' days'))
        Ani(k+1) = getframe;
    end

    if t(k+1)>=200-deltat/2 & t(k+1)<=200+deltat/2
        u200 = u;
    elseif t(k+1)>=600-deltat
        u600 = u;
    end

    difference = 1; %om de while-lus opnieuw te starten
end

%% Extra visualisatie(s)
figure(3);
% f.Position(3:4) = [1.5*560,420];
sgtitle(strcat('GS in xz-plane: y = 7 m'))
tijden = [200 600]
u_plot = {u200,u600};
id=find(y>(7-deltay/2) & y<(7+deltay/2))+1;
for p = 1:length(tijden)
    subplot(2,1,p)
    u = u_plot{p};
    u_xz = squeeze(u(id(1),2:end-1,2:end-1));
    contourf(x,z,u_xz')
    xlabel('x')
    ylabel('z')
    set (gca, 'ydir', 'reverse')
    title(strcat('t =',num2str(tijden(p)),'days'))
    colorbar
end
