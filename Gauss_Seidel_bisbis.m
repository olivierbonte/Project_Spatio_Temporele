%% Figuur 9 via Gauss-Seidel bis
clear all
close all
tic 
%Deze werkt!!!!!!
%% Initialisaties
v = 0.1; %m/d
alpha_L = 0.05; %m
D_L = alpha_L*v;
Pe = 1; %Peclet nummer
deltax = Pe*alpha_L; 

Cr = 1; %Courant nummer
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

u = zeros(N+4,M+4,K+1); %dus een 3D matrix
%initiele condities
idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
%beter idee! om rekening te houden met de imaginaire knopen
posx = find(idx) + 1; %omdat indexering start op -1
posy = find(idy) +1; %idem
C0 = 1; %bv 1 kg/m^3 dus
u(posy,posx,1) = C0;

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

%unew = zeros(N+4,M+4,K+1);
%unew = u;
u_flat = u(:,:,2); %laat deze dus op de initiële condities starten met z'n iteraties!
unew_flat = u_flat;
figure()
for k = 1:K %tot tijdstip K (K+1) maar indexeren op k+1
    iter = 0;
    while difference > phi
        for j = 2:N+3 % van j = 0 tot j = N+1
            for i = 2:M+3 % van i = 0 tot i = M+1
                unew_flat(j,i) = -B2/B1*u_flat(j,i+1)-B3/B1*unew_flat(j,i-1)-...
                    B4/B1*u_flat(j+1,i) - B5/B1*unew_flat(j-1,i) + R/(deltat*B1)*u(j,i,k);
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
    k
    u(:,:,k+1) = u_flat;
    %contourf(x,y,u(2:end-1,2:end-1,k+1));
    %xlabel('x')
    %ylabel('y')
    %zlabel('u')
    %title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    %colorbar
    %Ani(k) = getframe;
    difference = 1; %zodat hij terug vertrekt
end
%% visualisaties
%videostijl
figure()
for k = 1:10:K %Dus niet alle tijdstappen visualiseren
    contourf(x,y,u(2:end-1,2:end-1,k+1));
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title(strcat('tijd =  ',num2str(t(k+1))),'dagen')
    colorbar
    Ani(k) = getframe;
end

figure()
contourf(x,y,u(2:end-1,2:end-1,2))

figure()
contourf(x,y,u(2:end-1,2:end-1,3))

tijden = [0,200,600];
for i = 1:length(tijden)
    tijdid = t == tijden(i);
    
    figure()
    contourf(x,y,u(2:end-1,2:end-1,tijdid))
    xlabel('x')
    ylabel('y')
    title(strcat('t =',num2str(t(tijdid))))
end
toc
