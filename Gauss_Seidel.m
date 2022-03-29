%% Figuur 9 via Gauss-Seidel
clear all
close all
tic 
% FOUT: convergentie per tijdsitp is nodig!! 
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

%% Gauss Seidel iteratie
difference = 1;
phi = 1e-3; %convergentiecriterium
iter = 0;
%preallocaties
R = 1; %dus geen retardatie hier
B1 = R/deltat + 2*D_L/deltax^2 + 2*D_T/deltay^2;
B2 = -D_L/deltax^2 + v/(2*deltax);
B3 = -D_L/deltax^2 - v/(2*deltax);
B4 = -D_T/deltay^2;
B5 = B4;

unew = zeros(N+4,M+4,K+1);
unew(:,:,1) = u(:,:,1);

while difference > phi
iter = iter +1
    for k = 1:K %tot tijdstip K (K+1) maar indexeren op k+1
        for j = 2:N+3 % van j = 0 tot j = N+1
            for i = 2:M+3 % van i = 0 tot i = M+1
                unew(j,i,k+1) = -B2/B1*u(j,i+1,k+1)-B3/B1*unew(j,i-1,k+1)-B4/B1...
                    *u(j+1,i,k+1) - B5/B1*unew(j-1,i,k+1) + R/deltat*u(j,i,k);
                %bij index j-1 en i-1 dus unew omdat Gauss Seidel en niet
                %Jacobi!
            end
            unew(j,1,k+1) = unew(j,3,k+1); %Neumann links: u_-1 = u_1
            unew(j,M+4,k+1) = unew(j,M+2,k+1); %Neuman rehcts: u_M+2 = u_M
        end
        unew(1,:,k+1) = unew(3,:,k+1); %Neumann onderaan
        unew(N+4,:,k+1) = unew(N+2,:,k+1); %Neumann bovenaan
    end
    difference = max(max(max(abs(unew-u))))
    u = unew;
end
iter 

figure()
for i = 1:10
    contourf(x,y,u(2:end-1,2:end-1,i));
    xlabel('x')
    ylabel('y')
    zlabel('u')
    Ani(k) = getframe;
end


tijden = [0,200,600];
for i = 1:length(tijden)
    tijdid = t == tijden(i);
    figure()
    contourf(u(2:end-1,2:end-1,tijdid))
    xlabel('x')
    ylabel('y')
    title(strcat('t =',num2str(t(tijdid))))
end
toc
