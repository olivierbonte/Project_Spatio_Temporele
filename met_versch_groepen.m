%% Setup
clear all
clc
%indien niet vermeld: info uit paper 
%grid
dx = 50;    %km
x = 0:dx:6950;

dy = 50;     %km
y = 0:dy:12225;

dt = 0.1;      %dag+
t = 0:dt:20;

%parameters
K = 4.32; %km^2/dag    (50 m^2/s) 
mu = 0.7;   % %dag^-1
u = -100;   %km/dag   uit site: https://odnature.naturalsciences.be/marine-forecasting-centre/nl/maps/surface_sea_water_velocity/nos
v = -100;   %km/dag
k = 0;
phi = 0.1;
verschil = 10;

%intiele condities
N1 = ones(10, 10) *10000; %start van de larven populatie
N2 = ones(30,30) *10000; % """"
%voorbereiden van oplossing ruimte
N = zeros(length(y) + 2, length(x) + 2, length(t)); %2 imaginaire knopen
%dus N(y,x,t)
rijen_N = size(N,1);
kolommen_N = size(N,2);
a = [30, 50, 120, 50];
b = a+9;
c = [30,30, 120, 50];
d = c+9;

%de beginpopulaties plaatsen in het grote grid
for i = 1:4
    N(a(i):b(i),c(i):d(i),1) = N1;
end
N(210:239, 110:139, 1) = N2;

%
N(:,1,:) = 0;
N(:,kolommen_N,:)= 0; 
N(rijen_N,:,:)= 0;
N(1,:,:)= 0;





%initiele condities invoeren in oplossingsruimte
Nnew = N;
%loopt eerst volgens x-as dan volgens y-as en begint dan met het berekenen van
%oplossing voor volgend tijdstip


%% model runnen 
%de eerste dimensie van de matrix is de y-as
while verschil > phi
k = k+1
for i = 2:1:length(x)+1 %skipt de 2 imaginare knopen, x-as
    for g = 2:1:length(y)+1 %idem, y-as 
        for j = 1:1:length(t)-1 %t-as
            %1e afgeleide => voorwaarste differentie
            %2e afgeleide => centrale differentie
            Nnew(g,i,j+1) = dt*    (- u/dx * (N(g,i+1,j) - N(g,i,j)) - v/dy * (N(g+1,i,j) - N(g,i,j)) +...
                K * ( (N(g,i+1,j) - 2*N(g,i,j) + N(g, i-1,j))/dx^2 + (N(g+1,i,j) - 2*N(g,i,j) + N(g-1,i,j))/dy^2)) -mu*dt*N(g,i,j) + N(g,i,j);
        end
    end
end
verschil = max(max(max(abs(Nnew - N))));
N = Nnew;

end

%% Visualisatie 
figure(1)
%plotten voor bvb tijd = 10 dagen
tijd_num = 10;
opp = surf(y,x,N(2:rijen_N -1,2:kolommen_N - 1 ,tijd_num)');
title(['Concentratieprofiel op dag ', num2str(t(tijd_num)),'.'])
xlabel("y")
ylabel("x")
zlabel("concentratie aan larven")

figure(2)
tijd_num = 100;
opp = surf(y,x,N(2:rijen_N -1,2:kolommen_N - 1 ,tijd_num)');
title(['Concentratieprofiel op dag ', num2str(t(tijd_num)),'.'])
xlabel("y")
ylabel("x")
zlabel("concentratie aan larven")

figure(3)
x_p=zeros(1,length(x)+2);
y_p=zeros(1,length(y)+2); %Generate the plate
for i = 1:length(x)+2
    x_p(i) =(i-1)*dx;
end
for i = 1:length(y)+2
    y_p(i) =(i-1)*dy;
end
N_max = max(max(max(N)));
N_min = min(min(min(N)));


for j=1:201
    surf(y_p(1:end),x_p(1:end),N(:,:,j)')
    hold on
    title(sprintf('Concentration at days : %i days ',round(j*dt)))
    cb=colorbar;
    caxis([0 3]);
    view(90,-90);
    %xlim([0 6950+dx]); xlabel('Length');
    %ylim([0 12225+dy]); ylabel('Width');
    %zlim([0 3]); 
    zlabel('Concentratie');
    drawnow 
    hold off
end

Z = sum(sum(N));
a = zeros(1,size(N,3));
for i= 1:size(N,3)
    a(i) = Z(1,1,i);
end

figure(4)
plot(t, a)
ylabel("totale concentratie aan larven")
xlabel("tijd")
