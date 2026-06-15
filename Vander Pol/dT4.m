function [dT,dTx,dTy] = dT4(P,dy,dox,doy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function executes the finite differences method to solve numericaly 
% the Fokker-Planck equation.They were used approximations of 4th order for
% the derivatives (central-difference formula)
%% Input variables are described in FokkerPlanckEquation.m
%% Output variables: 
% - dT is a matrix that contains the values of the numerical solution of 
%   the Fokker-Planck equation
%% Application of the finite differences method
%  Coefficients corresponding to the approximation of the first derivative
% (of 4th order), of the Fokker-Planck equation solution
A1=[2 1 -1 -2; 4 1 1 4; 8 1 -1 -8; 16 1 1 16];                                                     
b1=[1; 0; 0; 0];
s1=A1\b1;
cf1=sum(s1);
%  Coefficients corresponding to the approximation of the second derivative
% (of 4th order), of the Fokker-Planck equation solution
A2=[3 2 1 -1 -2 -3; 9 4 1 1 4 9; 27 8 1 -1 -8 -27; 81 16 1 1 16 81; 243 32 1 -1 -32 -243; 729 64 1 1 64 729];
b2=[0; 1; 0; 0; 0; 0];
s2=A2\b2;
cf2=sum(s2);
% A will contain the coefficients of the system corresponding to the
% trivial solution of the stationary solution of the Fokker-Planck equation
% For this approximation it was used dx=dy 
A=zeros(length(dox)*length(doy));                                          % preallocated for efficiency                                      
% Coefficient for (x(i),y(j))                                           
for i=1:length(dox)
    for j=1:length(doy)
   A((i-1)*length(doy)+j,(i-1)*length(doy)+j)=-(cf2/(dy^2))+(cf1*doy(j)/dy)+(P(1)*dox(i)^2-P(2))-(cf1*((P(1)*dox(i)^2-P(2))*doy(j)+(P(4)^2)*dox(i))/dy);
    end
end
% Coefficient for(x(i),y(j+1))
for i=1:length(dox)
    for j=1:length(doy)-1
   A((i-1)*length(doy)+j,(i-1)*length(doy)+j+1)=(s2(3)/(dy^2))+(s1(2)*((P(1)*dox(i)^2-P(2))*doy(j)+(P(4)^2)*dox(i))/dy);
    end
end
% Coefficient for(x(i),y(j+2))
for i=1:length(dox)
    for j=1:length(doy)-2
   A((i-1)*length(doy)+j,(i-1)*length(doy)+j+2)=(s2(2)/(dy^2))+(s1(1)*((P(1)*dox(i)^2-P(2))*doy(j)+(P(4)^2)*dox(i))/dy);
    end
end
% Coefficient for(x(i),y(j+3))
for i=1:length(dox)
    for j=1:length(doy)-3
   A((i-1)*length(doy)+j,(i-1)*length(doy)+j+3)=(s2(1)/(dy^2));
    end
end
% Coefficient for(x(i),y(j-1))
for i=1:length(dox)
    for j=2:length(doy)
   A((i-1)*length(doy)+j,(i-1)*length(doy)+j-1)=(s2(4)/(dy^2))+(s1(3)*((P(1)*dox(i)^2-P(2))*doy(j)+(P(4)^2)*dox(i))/dy);
    end
end
% Coefficient for(x(i),y(j-2))
for i=1:length(dox)
    for j=3:length(doy)
   A((i-1)*length(doy)+j,(i-1)*length(doy)+j-2)=(s2(5)/(dy^2))+(s1(4)*((P(1)*dox(i)^2-P(2))*doy(j)+(P(4)^2)*dox(i))/dy);
    end
end
% Coefficient for(x(i),y(j-3))
for i=1:length(dox)
    for j=4:length(doy)
   A((i-1)*length(doy)+j,(i-1)*length(doy)+j-3)=(s2(6)/(dy^2));
    end
end
% Coefficient for(x(i-1),y(j))
for i=2:length(dox)
    for j=1:length(doy)
   A((i-1)*length(doy)+j,(i-2)*length(doy)+j)=-s1(3)*doy(j)/dy;
    end
end
% Coefficient for(x(i-2),y(j))
for i=3:length(dox)
    for j=1:length(doy)
   A((i-1)*length(doy)+j,(i-3)*length(doy)+j)=-s1(4)*doy(j)/dy;
    end
end
% Coefficient for(x(i+1),y(j))
for i=1:length(dox)-1
    for j=1:length(doy)
   A((i-1)*length(doy)+j,i*length(doy)+j)=-s1(2)*doy(j)/dy;
    end
end
% Coefficient for(x(i+2),y(j))
for i=1:length(dox)-2
    for j=1:length(doy)
   A((i-1)*length(doy)+j,(i+1)*length(doy)+j)=-s1(1)*doy(j)/dy;
    end
end
%% Coefficients of non-trivial linear system corresponding to the stationa-
%  ry solution of the Fokker-Planck equation
od=((1/2)*(length(dox)+1)-1)*length(doy)+((1/2)*(length(doy)+1));          % position of the central column of A
Ace=[A(:,1:(od-1)),A(:,(od+1):length(doy)*length(dox))];                   % Ace is the matrix A without its central column
Acfe=[Ace(1:(od-1),:);Ace((od+1):end,:)];                                  % Acfe is the matrix A without its central row
ce=A(:,od);                                                                % central column of A
Aoe=-[ce(1:(od-1));ce((od+1):end)];                                        % negative values of the central column of A,without the value corresponding to the central row of A
Sno=Acfe\Aoe;                                                              % solution of the linear system AcfeX=Aoe 
Snoc=[Sno(1:(od-1));1;Sno((od):end)];                                      % 
ISnoc=sum(Snoc*(dy^2));                                                    % normalization coefficient for estimation of the density (solution of the Fokker-Planck equation)
dT=zeros(length(doy),length(dox));
for i=1:length(dox)*length(doy)
    dT(i)=Snoc(i)/ISnoc;                                                   % numerical solution Fokker-Planck equation(approximation of invariant density)
end
% For this approximation it was used dx=dy 
dTx=sum(dy*dT);                                                            % marginal density with respect to x
dTy=sum(dy*dT');                                                           % marginal density with respect to y
end 

 