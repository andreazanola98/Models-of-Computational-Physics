%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    University of Padua, Department of Physics and Astronomy    %
%          Course of Computational Methods for Physics           %
%                  Lab session n°6 29/11/2019                    %
%                    Coded by: Andrea Zanola                     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
tic

%---Parameters of the System---
Lx=10;                  %x-Spatial Length
Nx=100;                 %n° x-points
hx=Lx/(Nx-1);           %Delta-x
x=0:hx:(Nx-1)*hx;       %Vector of possible x positions
Ly=10;
Ny=100;
hy=Ly/(Ny-1);
y=0:hy:(Ny-1)*hy;       

%---Charge Fixed Distribution---
px=4;                  %Down-Left x-coordinate position
py=4;                  %Down-Left y-"           "
Lc=1;                  %Length of the square lattice
xi=px*Nx/Lx +1;        %Initial x-index of the square lattice
xf=(px+Lc)*Nx/Lx +1;   %Final   x-index of the square lattice
yi=py*Ny/Ly +1;        %Initial y-index of the square lattice
yf=(py+Lc)*Ny/Ly +1;   %Final   y-index of the square lattice

rho0=100;   
F=zeros(Nx,Ny);                %F=rho/eps0, with eps0=1
F(yi:yi+2, xi-10:xi) = +rho0;  %Down Left  Corner
F(yf:yf+10, xi-2:xi) = -rho0;  %Up   Left  Corner
F(yf-2:yf, xf:xf+10) = +rho0;  %Up   Right Corner
F(yi-10:yi, xf:xf+2) = -rho0;  %Down Left  Corner


%---Gauss_Seidel's Method---
C = [1/hx^2 -2*(1/hx^2 +1/hy^2) 1/hy^2]; %Derivative operator
M_tilde = -[C(3) C(1) C(1) C(3)]'/C(2);

f_tilde=zeros(Nx*Ny,1,'single'); %Matrix to Vector rapresentation
for j=1:Ny
    for i=1:Nx
        l=(j-1)*Nx + i;
        f_tilde(l)=F(i,j)/C(2);
    end
end

PHI=f_tilde;
for z=1:10                  %z: number of iterations
    for j=1:Ny              %j: row index
        for i=1:Nx          %i: column index
            
            l=(j-1)*Nx+i;
            if i~=1 && i~=Nx && j~=1 && j~=Ny                       %If true the point is not on the border
                if z~=1                                             %If true then it's not the first iteration
                    phi_t = [PHI(l-Nx) PHI(l-1) PHI(l+1) PHI(l+Nx)];
                else
                    phi_t = [f_tilde(l-Nx) f_tilde(l-1) f_tilde(l+1) f_tilde(l+Nx)]; 
                end
                PHI(l) = phi_t*M_tilde + f_tilde(l);
            end
        end
    end 
end
 
phi=zeros(Nx, Ny);          %Vector to Matrix rapresentation
for j=1:Ny
    for i=1:Nx
        l=(j-1)*Nx+i;
        phi(i,j)=PHI(l);
    end
end

%---Electric Field Evaluation---
Ex=zeros(Nx, Ny);
Ey=zeros(Nx, Ny);
for j=1:Ny 
    for i=1:Nx
        if i~=1 && i~=Nx && j~=1 && j~=Ny   
            Ex(j,i)=0.5*(phi(j, i-1)-phi(j, i+1))/hx;
            Ey(j,i)=0.5*(phi(j-1, i)-phi(j+1, i))/hy;
        end
    end
end

toc

%---Plot Evaluation---
[X,Y]=meshgrid(x,y);

figure(1)
subplot(1,2,1)
surf(X, Y, F), title('Charge spatial distribution'),...
xlabel('x'), ylabel('y'), zlabel('\rho_0'), colorbar;
subplot(1,2,2)
surf(X, Y, phi), title('3D Potential'), xlabel('x'), ylabel('y'),...
zlabel('Potential V'), colorbar, grid minor;

figure(2)
quiver(X, Y, Ex, Ey,'LineWidth',0.75), hold on, contour(X, Y, F),...
title('Electric Field E'), xlabel('x'), ylabel('y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
