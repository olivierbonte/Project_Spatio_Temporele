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
deltay = deltax;
deltaz = deltax;
Cr=0.5;
deltat = Cr*deltax/v_x;

x = 0:deltax:120; 
y = 0:deltay:14; 
z = 2:deltaz:20; %18m dikke aquifer startend op 2m diepte (GWT)
t = 0:deltat:600; 
M = length(x)-2; 
N = length(y)-2;
O = length(z)-2;
K = length(t)-1;

u = zeros(N+4,M+4,O+4); %3D matrix, die telkens overschreven zal worden

% Initiële conditie
idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
idz = z <= 8; %in de bovenste 6m
posx = find(idx) + 1; %rekening houdend met imaginaire knopen
posy = find(idy) + 1;
posz = find(idz) + 1;
% Balk met C0=1 als initiële conditie?
C0 = 1; %bv 1 kg/m^3 dus
u(posy,posx,posz) = C0;

%Of bijv. lineaire afname tussen 2 en 8m diep 
for g=posz
    u(posy,posx,g)=1-1/6*(z(g-1)-2);
end

%% Plot initiele conditie
dieptes = [2,5,8];
for p=1:length(dieptes)
    id=find(z==dieptes(p))+1;
    figure(p)
    contourf(x,y,u(2:end-1,2:end-1,id))
    xlabel('x')
    ylabel('y')
    colorbar
end
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

u_k = u; %3D matrix van de vorige tijdstap (k) initialiseren
unew = u; %3D matrices van de huidige tijdstap (k+1) initialiseren

figure(4)
hold on
for k = 1:K %tot tijdstip K+1 door te indexeren op k+1
    iter = 0;
    while difference > phi
        for g = 2:O+3 %van g=0 tot g=O+1 
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
    if rem(t(k+1),5)==0
        sgtitle(strcat('GS in xy-plane, t =  ',num2str(t(k+1)),' days'));

        subplot(3,1,1)
        id1=find(z==dieptes(1))+1;
        contourf(x,y,u(2:end-1,2:end-1,id1));
        xlabel('x')
        ylabel('y')
        zlabel('C')
        title('z=2')
        colorbar

        subplot(3,1,2)
        id2=find(z==dieptes(2))+1;
        contourf(x,y,u(2:end-1,2:end-1,id2));
        xlabel('x')
        ylabel('y')
        zlabel('C')
        title('z=5')
        colorbar

        subplot(3,1,3)
        id3=find(z==dieptes(3))+1;
        contourf(x,y,u(2:end-1,2:end-1,id3));
        xlabel('x')
        ylabel('y')
        zlabel('C')
        title('z=8')
        colorbar

     Ani(k+1) = getframe;
    end

    if t(k+1) == 200
        u200 = u;
    elseif t(k+1) == 600
        u600 = u;
    end

    difference = 1; %om de while-lus opnieuw te starten
end

%% Extra visualisatie(s)
% Op t=200
figure(5);
hold on
sgtitle(strcat('GS in xy-plane, t = 200 days'));
for p=1:length(dieptes)
    id=find(z==dieptes(p))+1;
    subplot(3,1,p)
    contourf(x,y,u200(2:end-1,2:end-1,id))
    xlabel('x')
    ylabel('y')
    title(strcat('z =  ',num2str(dieptes(p))))
    colorbar
    caxis([0,0.25])
end
hold off

% Op t=600
figure(6);
hold on
sgtitle(strcat('GS in xy-plane, t = 600 days'));
for p=1:length(dieptes)
    id=find(z==dieptes(p))+1;
    subplot(3,1,p)
    contourf(x,y,u600(2:end-1,2:end-1,id))
    xlabel('x')
    ylabel('y')
    title(strcat('z =  ',num2str(dieptes(p))))
    colorbar
    caxis([0,0.1])
end
hold off
