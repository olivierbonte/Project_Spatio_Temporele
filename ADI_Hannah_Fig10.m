close all
clear all

%% Initialisations
% analoog Fig 10
v = 0.1; %m/d
alpha_L = 0.5; %m aanpssing tov Figuur 9, specifiek voor Figuur 10
D_L = alpha_L*v;
Pe = 0.5; %Peclet nummer naar 0.5
deltax = Pe*alpha_L; %Pe*alpha_L = Pe*D_L/v

Cr = 0.5;%Courant nummer 0.5
deltat = Cr*deltax/v; %in dagen. Hierbij is schema NIET stabiel
deltat = deltat/2; %omdat we in ADI per keer impliciet met een volledige deltat 
%vooruitgaan => nood aan delen door 2 voor stabiliteit! 
%experiment: deltatkleiner maken om te zien of golf trager gaat
%deltat = deltat/16;

alpha_T = 0.05; %m
D_T = alpha_T*v;

%grid
x = 0:deltax:120; %x tot 120 m
%bij gebrek aan info: neem hier deltay = deltax
deltay = deltax;
y = 0:deltay:14; %y tot 14m
t = 0:2*deltat:800; %Dus 2 deltats per volledige ADI doorlopen
M = length(x)-2; %aftrekken van interne knopen
N = length(y)-2;
%K = length(t)-1;
R = 1;

%% Coefficients Gaussian elimination
% Coefficienten voor impliciete x-afgeleide
alpha1=zeros(M+2,1); %alpha0 alpha1 ... alphaM+1
alpha1(1:end-1)=-D_L/deltax^2-v/(2*deltax);
alpha1(end)=-2*D_L/deltax^2; %alpha_tilde=gamma+alpha
beta1=R/deltat+D_L*2/deltax^2;
gamma1=zeros(M+2,1); %gamma0 gamma1 ... gammaM+1
gamma1(2:end)=-D_L/deltax^2+v/(2*deltax);
gamma1(1)=-2*D_L/deltax^2; %gamma_tilde=gamma+alpha

KappaCoeff1=zeros(M+2,1); %kappa-1 kappa0 ... kappaM
Y1=zeros(M+3,1); %Y-1 Y0 ... YM YM+1
RL1=zeros(M+2,1); %RL0 RL1 ... RLM+1
for i=2:M+2
    KappaCoeff1(i)=gamma1(i-1)/(beta1-alpha1(i-1)*KappaCoeff1(i-1));
end

% Coefficienten voor impliciete y-afgeleide
alpha2=zeros(N+2,1); %alpha0 alpha1 ... alphaN+1
alpha2(1:end-1)=-D_T/deltay^2;
alpha2(end)=-2*D_T/deltay^2; %alpha_tilde=gamma+alpha
beta2=R/deltat+D_T*2/deltay^2;
gamma2=zeros(N+2,1); %gamma0 gamma1 ..gammaN+1
gamma2(2:end)=-D_T/deltay^2;
gamma2(1)=-2*D_T/deltay^2; %gamma_tilde=gamma+alpha

KappaCoeff2=zeros(N+2,1); %kappa-1 kappa0 ... kappaN
Y2=zeros(N+3,1); %Y-1 Y0 ... YN+1
RL2=zeros(N+2,1); %RL0 RL1 ... RLN+1
for j=2:N+2
    KappaCoeff2(j)=gamma2(j-1)/(beta2-alpha2(j-1)*KappaCoeff2(j-1));
end

C=zeros(N+2,M+2);

%% Initial condition
%for i=find(x==10):find(x==14.600000000000001)
%    for j=find(y==4.2):find(y==9.8)
%        C(j,i)=1; %NIET GEGEVEN!!!
%    end
%end
%altneratieve implementatie, analoog aan Fig10_GS script
idx = x >= 10 & x <= 14.6;
idy = y >= 4.2 & y <= 9.8;
%beter idee! om rekening te houden met de imaginaire knopen
posx = find(idx); 
posy = find(idy); 
C0 = 1; %bv 1 kg/m^3 dus
C(posy,posx) = C0;

%% Plot of the initial condition
step=0;
surf(x,y,C,'linestyle','none');
xlabel('x')
ylabel('y')
zlabel('c')
title('Project 1')
Ani(1) = getframe;


%% ADI
figure()
for tijd = t
    step=step+1;
    for j=1:N+2
        for i=1:M+2
            if j==1
                RL1(i)=(R/deltat-D_T*2/deltay^2)*C(j,i)+2*D_T/deltay^2*C(j+1,i); %niet meer steunend op C(j-1,i)
            elseif j==N+2
                RL1(i)=2*D_T/deltay^2*C(j-1,i)+(R/deltat-D_T*2/deltay^2)*C(j,i); %niet meer steunend op C(j+1,i)
            else
                RL1(i)=D_T/deltay^2*C(j-1,i)+(R/deltat-D_T*2/deltay^2)*C(j,i)+D_T/deltay^2*C(j+1,i);
            end
            Y1(i+1)=(RL1(i)-alpha1(i)*Y1(i))/(beta1-alpha1(i)*KappaCoeff1(i));
        end
        C(j,M+2)=Y1(M+3);
        for i=M+1:-1:1
            C(j,i)=Y1(i+1)-KappaCoeff1(i+1)*C(j,i+1);
        end
    end

    for i=1:M+2
        for j=1:N+2
            %aanpassing in komma's geen effect logisherwijs (termen
            %annuleren elkaar)
            if i==1
                RL2(j)=(R/deltat-D_L*2/deltax^2)*C(j,i)+2*D_L/deltax^2*C(j,i+1);
                %RL2(j)=(D_L/deltax^2+v/(2*deltax))*C(j,i+1)+(R/deltat-D_L*2/deltax^2)*C(j,i)+(D_L/deltax^2-v/(2*deltax))*C(j,i+1);
            elseif i==M+2
                RL2(j)=2*D_L/deltax^2*C(j,i-1)+(R/deltat-D_L*2/deltax^2)*C(j,i);
                %RL2(j)=(D_L/deltax^2+v/(2*deltax))*C(j,i-1)+(R/deltat-D_L*2/deltax^2)*C(j,i)+(D_L/deltax^2-v/(2*deltax))*C(j,i-1);
            else
                RL2(j)=(D_L/deltax^2+v/(2*deltax))*C(j,i-1)+(R/deltat-D_L*2/deltax^2)*C(j,i)+(D_L/deltax^2-v/(2*deltax))*C(j,i+1);
            end
            Y2(j+1)=(RL2(j)-alpha2(j)*Y2(j))/(beta2-alpha2(j)*KappaCoeff2(j));
        end
        C(N+2,i)=Y2(N+3);
        for j=N+1:-1:1
            C(j,i)=Y2(j+1)-KappaCoeff2(j+1)*C(j+1,i);
        end
    end
    %surf(x,y,C,'linestyle','none');
    if rem(tijd,5) == 0
        contourf(x,y,C)
        colorbar
        xlabel('x')
        ylabel('y')
        zlabel('C')
        title(strcat('ADI: time = ',num2str(tijd), 'days'))
        Ani(step+1) = getframe;
    end
    %voor fig 10 nabootsen
    if tijd == 200
        C200 = C;
    elseif tijd == 600
        C600 = C;
    end
end

%% Fig 10 nabootsen
tijden = [200,600];
f = figure();
f.Position(3:4) = [1.5*560,420];
Cs = {C200,C600};
for i = 1:length(tijden)
    subplot(2,1,i)
    C = Cs{i};
    contourf(x,y,C)
    xlabel('x')
    ylabel('y')
    title(strcat('t =',num2str(tijden(i)),'days'))
    colorbar
end
sgtitle(strcat('ADI, \Deltat = ',num2str(deltat),' days'))
exportgraphics(gcf,'Figuur_10_ADI.png','Resolution',900)
toc


%% Play the animation

%movie(Ani,0)

video = VideoWriter('With_diff_movie.mp4');
video.FrameRate = 60;
open(video);
writeVideo(video, Ani);
close(video);