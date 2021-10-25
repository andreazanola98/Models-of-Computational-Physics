%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for spin flipping and evalutation of a new system     %
% spin configuration.                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H, M, sig] = CONF_SPIN(NL, h, J, b, sig)

%---Sping Flipping---
if b==1                 
     xx=randi([1 NL],1,1);
     xy=randi([1 NL],1,1);
     
     sig(xx, xy) = -sig(xx, xy);
     
     %---Updating BC
     if xx==1 && xy~=1
        sig(NL+1, xy) = -sig(NL+1, xy);
     end
     if xy==1 && xx~=1
        sig(xx, NL+1) = -sig(xx, NL+1);
     end
else
    %At the beginning b=0, because the matrix isn't built yet.
    sig=int8(2*randi([0 1], NL+1, NL+1)-1);
end

%---Boundary Conditions---
sig(1:NL, 1) = sig(1:NL, NL+1); %1° column is a copy of the last one
sig(1, 1:NL) = sig(NL+1, 1:NL); %1° raw    "     "           "
sig(NL+1, NL+1)=0;

%---Physical Quantities Calculation---
M=sum(sum(sig(1:NL, 1:NL)));  %Total magnetization
Hsig=h*M;                     %Total magnetic field

HJ=0;                         %Coupling energy inizialization
for i=1:NL
    for j=1:NL
        if sig(i,j)*sig(i,j+1)==1 %Horizontal interactions
            HJ=HJ-J;
        else
            HJ=HJ+J;
        end
        
        if sig(j,i)*sig(j+1,i)==1 %Vertical interactions
            HJ=HJ-J;
        else
            HJ=HJ+J;
        end
        
    end
end

H=-Hsig+HJ; %Total Energy

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
