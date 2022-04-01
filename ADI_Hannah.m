close all
clear all

%% Initialisations
 deltax=0.05;
 deltay=0.05; 
 deltat=0.1; 
 x=0:deltax:120;
 y=0:deltay:14;
 M=length(x)-2;
 N=length(y)-2;

 v=0.1; %m/d
 alpha_L=0.05; %m
 D_L=alpha_L*v; %m2/s
 alpha_T=0.005;
 D_T=alpha_T*v;
 R=1;

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
 for i=find(x==10):find(x==14.600000000000001)
     for j=find(y==4.2):find(y==9.8)
         C(j,i)=1; %NIET GEGEVEN!!!
     end
 end

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
for t=0:deltat:800
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
             if i==1
                 RL2(j)=(R/deltat-D_L*2/deltax^2)*C(j,i)+2*D_L/deltax^2*C(j,i+1);
             elseif i==M+2
                 RL2(j)=2*D_L/deltax^2*C(j,i-1)+(R/deltat-D_L*2/deltax^2)*C(j,i);
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
    if rem(t,5) == 0
        contourf(x,y,C)
        colorbar
        xlabel('x')
        ylabel('y')
        zlabel('C')
        title(strcat('ADI: time = ',num2str(t)))
        Ani(step+1) = getframe;
    end
end
%% Play the animation

%movie(Ani,0)

video = VideoWriter('With_diff_movie.mp4');
video.FrameRate = 60;
open(video);
writeVideo(video, Ani);
close(video);