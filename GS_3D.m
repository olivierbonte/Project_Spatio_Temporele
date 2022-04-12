%% Parameters van figuur 10, nu in 3D
%clear all
close all
tic 
%Deze werkt!!!!!!
%% Initialisaties
v = 0.1; %m/d
alpha_L = 0.5; %m aanpssing tov Figuur 9
D_L = alpha_L*v;
Pe = 0.5; %Peclet nummer, moet kleiner dan 1
deltax = Pe*alpha_L; 

Cr = 0.5; %Courant nummer
deltat = Cr*deltax/v; %in dagen

alpha_T = 0.05; %m
D_T = alpha_T*v; 

%z richting
alpha_V = 0.05; %V van verticaal dus
D_V = alpha_V*v;

%grid
x = 0:deltax:120; %x tot 120 m
%bij gebrek aan info: neem hier deltay = deltax
deltay = deltax;
deltaz = deltax; %analoog dus
y = 0:deltay:14; %y tot 14m
z = 2:deltaz:20; 
%idee van GWT op 2 m
t = 0:deltat:600; %tot dag 600 in de figuur
M = length(x)-2; %aftrekken van interne knopen
N = length(y)-2;
Z = length(z)-2;
K = length(t)-1;

u = zeros(N+4,M+4,Z+4); %dus een 3D matrix!!
%initiele condities
idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
idz = z >= 2 & z <= 8; 
posx = find(idx) + 1; %omdat indexering start op -1
posy = find(idy) +1; %idem
posz = find(idz); %hier niet selecteren uit de matrix
concentraties = linspace(1,0,length(posz));
for i = 1:length(posz)
    C0 = concentraties(i);
    u(posy,posx,posz(i)) = C0;
end

%plot de initiele conditie in xz vlak
ydoorsnede = (4.2 + 9.8)/2;
idy_doorsnede = y == ydoorsnede;
posy_doorsnede = find(idy_doorsnede) +1;
figure()
uxz = squeeze(u(posy_doorsnede,2:end-1,2:end-1));
contourf(x,z,uxz')
%gebruik squeeze om 2D matrix te krijgen!!!
xlabel('x [m]')
ylabel('z [m]')
set ( gca, 'ydir', 'reverse' )
colorbar
%title(strcat('Initial conditions in the xz-plane: y = ',num2str(ydoorsnede)',' m'))
title(strcat('y = ',num2str(ydoorsnede)',' m'))
exportgraphics(gcf,'Init_cond_3D.png','Resolution',900)
%% Gauss Seidel iteratie
difference = 1;
phi = 1e-4; %convergentiecriterium
iter = 0;
%preallocaties
R = 1; %dus geen retardatie hier
B1 = R/deltat + 2*D_L/deltax^2 + 2*D_T/deltay^2 + 2*D_V/deltaz^2;
B2 = -D_L/deltax^2 + v/(2*deltax);
B3 = -D_L/deltax^2 - v/(2*deltax);
B4 = -D_T/deltay^2;
B5 = B4;
B6 = -D_V/deltaz^2;
B7 = B6;

unew = u;
u_previous = u; %Dit is de vorige tijdstap altijd opgeslaan
figure()
for k = 1:K %indexeer op k+1
    iter = 0;
    while difference > phi
        for j = 2:N+3 % van j = 0 tot j = N+1
            for i = 2:M+3 % van i = 0 tot i = M+1
                for g = 2:Z+3 
                unew(j,i,g) = -B2/B1*u(j,i+1,g)-B3/B1*unew(j,i-1,g)-...
                    B4/B1*u(j+1,i,g) - B5/B1*unew(j-1,i,g) - ...
                    B6/B1*u(j,i,g+1) - B7/B1*unew(j,i,g-1) +...
                    R/(deltat*B1)*u_previous(j,i,g);
                end
                %ook in z richting 2 neumann randvoorwaarden!
                unew(j,i,Z+4) = unew(j,i,Z+2);
                unew(j,i,1) = unew(j,i,3);
            end
            unew(j,1,g) = unew(j,3,g); %Neumann links: u_-1 = u_1
            unew(j,M+4,g) = unew(j,M+2,g); %Neuman rehcts: u_M+2 = u_M
        end
        unew(1,:,g) = unew(3,:,g); %Neumann onderaan
        unew(N+4,:,g) = unew(N+2,:,g); %Neumann bovenaan
        difference = max(max(max(abs(unew-u))));
        u = unew;
        iter = iter +1;
    end
    u = unew;
    u_previous = unew; %dus bijhouden als vorige tijdstap altijd! 
    if rem(t(k+1),10) == 0 %dus niet alle tijdstappen visualiseren!
        k
        uxz = squeeze(u(posy_doorsnede,2:end-1,2:end-1));
        contourf(x,z,uxz')
        %gebruik squeeze om 2D matrix te krijgen!!!
        xlabel('x')
        ylabel('z')
        set ( gca, 'ydir', 'reverse' )
        colorbar
        title(strcat('GS in xz-plane: y = ',num2str(ydoorsnede)',' m, t = ',...
            num2str(t(k+1)),' days'))
        Ani(k+1) = getframe;
    end
    if t(k+1) == 200
        u200 = u;
    elseif t(k+1) == 600
        u600 = u;
    end
    difference = 1; %zodat hij terug vertrekt
end
%% visualisaties
%videostijl
%figure()
%{
for k = 0:1:K %Dus niet alle tijdstappen visualiseren
    contourf(x,y,u10(2:end-1,2:end-1,k+1));
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    colorbar
    Ani(k+1) = getframe;
end
%}

% xz vlak
tijden = [200,600];
f = figure();
f.Position(3:4) = [1.5*560,420];
us = {u200,u600};
for i = 1:length(tijden)
    tijdid = t == tijden(i);
    subplot(2,1,i)
    u = us{i};
    uxz = squeeze(u(posy_doorsnede,2:end-1,2:end-1));
    contourf(x,z,uxz')
    xlabel('x [m]')
    ylabel('z [m]')
    set ( gca, 'ydir', 'reverse' )
    title(strcat('t =',num2str(t(tijdid)),' days'))
    colorbar
    caxis([0,0.25])
end
sgtitle(strcat('xz-plane: y = ',num2str(ydoorsnede)',' m'))
exportgraphics(gcf,'Figuur_3D_xz_y7.png','Resolution',900)
toc

% xy vlak doorsnedes op t = 200
f = figure();
f.Position(3:4) = [1.5*560,420];
zdoorsnedes = [2,5,8];
for i = 1:length(zdoorsnedes)
    zdoorsnede = zdoorsnedes(i);
    idz_doorsnede = z == zdoorsnede;
    posz_doorsnede = find(idz_doorsnede) +1;
    uyx = squeeze(u200(2:end-1,2:end-1,posz_doorsnede));
    subplot(3,1,i)
    contourf(x,y,uyx)
    set ( gca, 'ydir', 'reverse' )
    xlabel('x [m]')
    ylabel('y [m]')
    colorbar
    caxis([0,0.25])
    title(strcat('z = ',num2str(zdoorsnede), ' m'))
end
sgtitle('xy-plane: t = 200 days')

% xy vlak doorsende op t = 600
f = figure();
f.Position(3:4) = [1.5*560,420];
zdoorsnedes = [2,5,8];
for i = 1:length(zdoorsnedes)
    zdoorsnede = zdoorsnedes(i);
    idz_doorsnede = z == zdoorsnede;
    posz_doorsnede = find(idz_doorsnede) +1;
    uyx = squeeze(u600(2:end-1,2:end-1,posz_doorsnede));
    subplot(3,1,i)
    contourf(x,y,uyx)
    set ( gca, 'ydir', 'reverse' )
    colorbar
    xlabel('x [m]')
    ylabel('y [m]')
    caxis([0,0.1])
    title(strcat('z = ',num2str(zdoorsnede), ' m'))
end
sgtitle('xy-plane: t = 600 days')