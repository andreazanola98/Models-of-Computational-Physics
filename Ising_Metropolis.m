%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    University of Padua, Department of Physics and Astronomy    %
%          Course of Computational Methods for Physics           %
%                  Lab session n°7 13/12/2019                    %
%                    Coded by: Andrea Zanola                     % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
close all
tic

%---Parameters of the system---
NL=20;         %Dimensions of the squared lattice (NLxNL)
h=0.001;           %External magnetic field
J=1;           %Ferromagnetic coupling strength
kbT=1.8;         %kb*T value [Remember: kbTc= 2/ln(1+sqrt(2))= 2.27

%---Quantity of Interest---
N=10^6;         %Monte Carlo steps
Hv=zeros(1, N); %Energy (H=Hamiltonian)
Mv=zeros(1, N); %Magnetization x Site

D=10;                        %Number of evalutaions for fluctuations
MT=zeros(1, length(kbT));    %Magnetization for various temperatures
devMT=zeros(1, length(kbT)); %Standard Dev. "      "          " 
MTi=zeros(1, D);             %Magnetization for various evalutations
devMTi=zeros(1, D);          %Standard Dev. "      "          " 

%---Metropolis Algorithm---
sig=zeros(NL+1, 'int8');   %Configuration of the system (we add +1 for BC)
[H0, M0, Sig] = CONF_SPIN(NL, h, J, 0, sig); 

for k=1:N
    [H1, M1, sig] = CONF_SPIN(NL, h, J, 1, Sig); 
    Eflip=cast(H1-H0,'double');
    r=rand/exp(-Eflip/kbT); 
 
    if Eflip<=0 || r<=1
        Sig=sig;
        M0=M1;
        H0=H1;
        Hv(k)=H1;
        Mv(k)=M1;
    else
        Hv(k)=H0;
        Mv(k)=M0;
    end
end

Mv=Mv/NL^2;
Hv=Hv/NL^2;

toc

%---Plot Visualization---
figure(1)
subplot(1,2,1)
plot(1:1:N, Mv,'r'), title('Magnetization x Site'), 
xlabel('Monte Carlo Steps'), ylabel('m'), axis([1 N -1 1]), grid on;
subplot(1,2,2)
plot(1:1:N, Hv,'r'), title('Energy x Site'), 
xlabel('Monte Carlo Steps'), ylabel('H'), grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




