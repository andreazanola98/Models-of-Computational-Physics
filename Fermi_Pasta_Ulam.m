%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    University of Padua, Department of Physics and Astronomy    %
%          Course of Computational Methods for Physics           %
%                  Lab session n°4 15/11/2019                    %
%                    Coded by: Andrea Zanola                     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
tic

%---Parameters of the system---
N=31;               %Number of masses
ni=1:N; 
n=1:N+2;            %N masses + 2 fixed nodes (walls)

m=ones(1, N+2);     %Masses
k=ones(1, N+1);     %Elasticity constants
a=0.25;             %Non linear parameter Alfa
J=1;                %Foundamental node
K=pi/(N+1);

%---Parameters of the simulation---
dt=0.1;      %[s]
tmax=10000;  %[s]
t=0:dt:tmax;
lt=length(t);

%---Initial Conditions---
x=zeros(N+2, lt);       
x(2:N+1, 1)=sin(K*ni*J); %Initial elongations
v=zeros(N+2, lt);        %Initial velocity
f=zeros(1, N);           %Initial forces

%---Implementation Velocity Verlet---
for j=1:lt-1 %j: time index
    
    for i=2:N+1 %i: atom index
        f(i-1) = (k(i)*x(i+1,j) - x(i,j)*(k(i)+k(i-1)) + k(i-1)*x(i-1,j))/m(i)...
                 + a*(-x(i-1,j)^2 + 2*x(i-1,j)*x(i,j) - 2*x(i,j)*x(i+1,j) + x(i+1,j)^2 );
    end 
 
    for i=2:N+1
        x(i,j+1) = x(i,j) + v(i,j)*dt + 0.5*f(i-1)*dt^2/m(i);
    end
 
    for i=2:N+1 
         f(i-1)= f(i-1)+(k(i)*x(i+1,j+1) - x(i,j+1)*(k(i)+k(i-1)) + k(i-1)*x(i-1,j+1))/m(i)...
                + a*(-x(i-1,j+1)^2 + 2*x(i-1,j+1)*x(i,j+1) - 2*x(i,j+1)*x(i+1,j+1) + x(i+1,j+1)^2 ); 
    end
 
    for i=2:N+1
        v(i,j+1) = v(i,j) + 0.5*f(i-1)*dt/m(i);
    end
 
end

%---Energy Calculations---
Ap=5;             %Select the number of nodes of interest
E=zeros(Ap, lt);  %Energy vector
C2=zeros(1, lt);  %Cn(t) is the contribution to the motion at node n
dC2=zeros(1, lt); %dCn(t) is the temporal derivative of Cn(t)

for z=1:Ap %z: mode index
    C2  = (sin(K*ni*z)*x(2:N+1,:)).^2;
    dC2 = (sin(K*ni*z)*v(2:N+1,:)).^2;

    for i=1:N  %sum over the number of nodes
        E(z,:) = E(z,:) + 0.5*C2*(k(i+1)*sin(K*z*(i+1)) - k(i)*sin(K*z*i))^2 + 0.5*m(i)*dC2*sin(K*z*i)^2;
    end
end

Et=sum(E); %Total Energy
toc

%---Plot Visualization---
figure(1)
for i=1:Ap
   plot(t,E(i,:)), hold on;
end
plot(t, Et,'k')
title({'Problema di Fermi-Pasta-Ulam','Energia Modi Normali'}),...
xlabel('Time'), ylabel('E'); 
legend('Node 1','Node 2','Node 3','Node 4','Node 5','E_{TOT}');
ax = gca; ax.FontSize = 20; ax.LineWidth=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

