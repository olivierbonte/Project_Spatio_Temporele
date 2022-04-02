%% Parameters van figuur 10, nu in 3D
%clear all
%close all
tic 
%Deze werkt!!!!!!
%% Initialisaties
v = 0.1; %m/d
alpha_L = 0.5; %m aanpssing tov Figuur 9
D_L = alpha_L*v;
Pe = 1; %Peclet nummer
deltax = Pe*alpha_L; 

Cr = 1; %Courant nummer
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
xlabel('x')
ylabel('z')
set ( gca, 'ydir', 'reverse' )
colorbar
title(strcat('Initial conditions in the xz-plane: y = ',num2str(ydoorsnede)',' m'))

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
%B6 = 
%B7 = 

%unew = zeros(N+4,M+4,K+1);
%unew = u;
u_flat = u10(:,:,2); %laat deze dus op de initiële condities starten met z'n iteraties!
unew_flat = u_flat;
figure()
for k = 1:K %tot tijdstip K (K+1) maar indexeren op k+1
    iter = 0;
    while difference > phi
        for j = 2:N+3 % van j = 0 tot j = N+1
            for i = 2:M+3 % van i = 0 tot i = M+1
                unew_flat(j,i) = -B2/B1*u_flat(j,i+1)-B3/B1*unew_flat(j,i-1)-...
                    B4/B1*u_flat(j+1,i) - B5/B1*unew_flat(j-1,i) + R/(deltat*B1)*u10(j,i,k);
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
    iter;
    k;
    u10(:,:,k+1) = u_flat;
    %contourf(x,y,u(2:end-1,2:end-1,k+1));
    %xlabel('x')
    %ylabel('y')
    %zlabel('u')
    %title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    %colorbar
    %Ani(k) = getframe;
    if rem(k,5) == 0
       k 
    end
    difference = 1; %zodat hij terug vertrekt
end
%% visualisaties
%videostijl
figure()
for k = 0:1:K %Dus niet alle tijdstappen visualiseren
    contourf(x,y,u10(2:end-1,2:end-1,k+1));
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    colorbar
    Ani(k+1) = getframe;
end

tijden = [200,600];
f = figure();
f.Position(3:4) = [1.5*560,420];
for i = 1:length(tijden)
    tijdid = t == tijden(i);
    subplot(2,1,i)
    contourf(x,y,u10(2:end-1,2:end-1,tijdid))
    xlabel('x')
    ylabel('y')
    title(strcat('t =',num2str(t(tijdid))),'days')
    colorbar
end
exportgraphics(gcf,'Figuur_10.png','Resolution',900)
toc