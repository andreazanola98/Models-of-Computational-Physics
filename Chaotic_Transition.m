%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    University of Padua, Department of Physics and Astronomy    %
%          Course of Computational Methods for Physics           %
%                    Final Project 30/01/2020                    %
%                    Coded by: Andrea Zanola                     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
tic 

%---Simulation Parameters---
dt=0.01;       %[s]
tmax=1000;     %[s]
t=0:dt:tmax;
N=length(t)-1;

t_reduced=100/dt;

%---Pendulum Parameters---
g=9.81;     %Gravity
l=9.81;     %Length of pendulum
m=1;        %Mass
k=g/l;

%---Dynamical Parameters---
theta0=0.2; %Initial angle [rad]
w0=0;       %Initial angular velocity [rad/s]
q=0.5;      %Damping
Fd=0.2;     %Amplitude of the external force
wd=2/3;     %Angular velocity of the external force

%---Equations of Motion (approx. sin(theta)~theta)---
Tap=2*pi*sqrt(1/k);               %Period
w=sqrt(k);                        %Intrinsic angualar velocity

fi0 = atan(-w0/(w*theta0));        %Initial phase
A = sqrt(theta0^2 + (w0/w)^2);     %Amplitude of oscillations
thetap = (180/pi)*A*cos(w*t+fi0);  %Trajectory
wp = -A*w*sin(w*t+fi0);            %Angular velocity

%---Runge Kutta 4 Integration---
Y=zeros(2,N);
Y(:,1)=[theta0; w0];
for i=1:N-1
    
  Y1 = Y(:,i);
  fY1= [Y1(2) ; -k*sin(Y1(1)) - q*Y1(2) + Fd*sin(wd*i*dt)]; 
  Y2 = Y1 + fY1*(dt/2);
  fY2= [Y2(2) ; -k*sin(Y2(1)) - q*Y1(2) + Fd*sin(wd*i*dt)]; 
  Y3 = Y1 + fY2*(dt/2);
  fY3= [Y3(2) ; -k*sin(Y3(1)) - q*Y1(2) + Fd*sin(wd*i*dt)]; 
  Y4 = Y1 + fY3*(dt);
  fY4= [Y4(2) ; -k*sin(Y4(1)) - q*Y1(2) + Fd*sin(wd*i*dt)];
  
  Y(:,i+1) = Y1 + (fY1 +2*fY2 +2*fY3 + fY4)*(dt/6);
  
  %---Keeping the angle inside 0-pi---
  if Y(1,i+1)>pi
    Y(1,i+1) = Y(1,i+1) -2*pi;
  elseif Y(1,i+1)<-pi
    Y(1,i+1) = Y(1,i+1) +2*pi;
  end
  
end
theta=Y(1,:)*(180/pi);

%---Cartesian Coordinates---
X=l*sin(Y(1,:));    
Z=l*(1-cos(Y(1,:)));

Xap=l*sin(thetap*(pi/180));
Zap=l*(1-cos(thetap*(pi/180)));

%---Energy Calculations---
E=0.5*m*(Y(2,:)*l).^2 +m*g*Z;
Eap=0.5*m*(wp*l).^2 +m*g*Zap;

%---Strange Attractor---
J=zeros(1,N);
for i=1:N
    for a=1:2000
        if abs(t(i) - 2*a*pi/wd) < dt/2
            J(i)=i;
            break
        end
    end
end
J(J==0)=[];

attractor=zeros(2,length(J));
for i=1:length(J)
    attractor(:,i)=Y(:,J(i));
end

%---Initial Condition Dependence Chaotic Motion---
delta_theta=0.001;
theta1=theta0 + delta_theta;
W=zeros(2,N);
W(:,1)=[theta1; w0];
for i=1:N-1
    
  Y1 = W(:,i);
  fY1= [Y1(2) ; -k*sin(Y1(1)) - q*Y1(2) + Fd*sin(wd*i*dt)]; 
  Y2 = Y1 + fY1*(dt/2);
  fY2= [Y2(2) ; -k*sin(Y2(1)) - q*Y1(2) + Fd*sin(wd*i*dt)]; 
  Y3 = Y1 + fY2*(dt/2);
  fY3= [Y3(2) ; -k*sin(Y3(1)) - q*Y1(2) + Fd*sin(wd*i*dt)]; 
  Y4 = Y1 + fY3*(dt);
  fY4= [Y4(2) ; -k*sin(Y4(1)) - q*Y1(2) + Fd*sin(wd*i*dt)];
  
  W(:,i+1)= Y1 + (fY1 +2*fY2 +2*fY3 + fY4)*(dt/6);
  
  %---Keeping the angle inside 0-pi---
  if W(1,i+1)>pi
    W(1,i+1) = W(1,i+1) -2*pi;
  elseif W(1,i+1)<-pi
    W(1,i+1) = W(1,i+1) +2*pi;
  end
  
end
psi=W(1,:)*(180/pi);
delta=abs(theta-psi);

toc

%---Visualization---
figure(1) %Strange Attractor
scatter(attractor(1,:)*(180/pi), attractor(2,:),'r.');
title('Pendulum Strange Attractor'), xlabel('\theta [°]'), ylabel('\omega [rad/s]'),...
grid minor, axis([-190 190 -2 1]);

figure(2) %Dynamics
subplot(2,2,1)
plot(t(1:t_reduced), theta(1,1:t_reduced),'r'), hold on, plot(t(1:t_reduced),thetap(1:t_reduced),'b');
title('Trajectory'), xlabel('t [s]'), ylabel('\theta [°]'), grid on, axis([0 max(t_reduced*dt) -180 180]);

subplot(2,2,2)
plot(t(1:t_reduced), Y(2,1:t_reduced),'r'),hold on, plot(t(1:t_reduced),wp(1:t_reduced),'b');
title('Angular Velocity'), xlabel('t [s]'), ylabel('\omega [rad/s]'), grid on;
 
subplot(2,2,3)
plot(theta(1,1:t_reduced), Y(2,1:t_reduced),'r'),hold on, plot(thetap,wp,'b');
title('Phase Portrait'), xlabel('\theta [°]'), ylabel('\omega [rad/s]'), grid minor;
legend('Sol. RK -4', 'Sol. sin(\theta)~\theta');

subplot(2,2,4)
plot(t(1:t_reduced), E(1:t_reduced),'r'), hold on, plot(t(1:t_reduced),Eap(1:t_reduced),'b'); 
title('Energy'), xlabel('t [s]'), ylabel('E [J]'), grid on

figure(3) %Different Initial Conditions
semilogy(t(1:2*t_reduced), delta(1:2*t_reduced), 'r.');
title('Orbits Deviation'), ylabel('\Delta\theta [°]'), xlabel('t [s]'), grid on;

figure(4) %Lissjous figure
comet(theta(1:t_reduced),thetap(1,1:t_reduced));
title('Lissajous figure'), xlabel('RK-4 Solution'),ylabel('sin(\theta) ~\theta Solution');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

