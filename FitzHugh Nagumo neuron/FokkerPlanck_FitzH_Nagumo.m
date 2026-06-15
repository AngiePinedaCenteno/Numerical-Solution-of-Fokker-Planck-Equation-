function FokkerPlanck_FitzH_Nagumo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a numerical solution of the Fokker-Planck equa-
% tion corresponding to the FitzHugh-Nagumo model in R^2:
%
% dX_t=Y_tdt
% dY_t=-sigma dW_t-(1/epsilon)*[(3X_t^2-1+epsilon)Y_t+ X_t^3+(gamma-1)X_t^2+ s + beta]dt
%
% where sigma=\tilde(sigma)/epsilon.
%% Output variables: 
% - dT contains the values of the numerical solution of the Fokker-Planck 
%   equation

%% Problem parameters: 
% - dx,dy:grid spacings for x and y in the discretized domain of density
% - dox: discretized domain in x direction for the density
% - doy: discretized domain in y direction for the density
% - PFH: vector that contains the parameters of the model defined above. 
%   Specifically, PFH=(epsilon,gamma,s,beta,\tilde(sigma))
%%
%%
tic
%% Assignation of the parameters values
PFH=[1.5,0.2,0.3,2,0.9];                                                     
P=zeros(1,6);
P(1)=3/PFH(1);P(2)=(1/PFH(1))-1; P(3)=(PFH(5)/PFH(1))^2;P(4)=PFH(2)-1; P(5)=PFH(3)+PFH(4);P(6)=1/PFH(1);

%% Discretization of the  domain of the density. Denominator of f must be a
%% pair number
dy=1/12;dx=dy;dox=(-3:dy:3);doy=(-4.5:dy:4.5);                                   

%% Remark: 
% To obtain the numerical solution of the Fokker-Planck equation it was 
% used  the finite differences method. The grid used to compute the numeri-
% cal solution  is a rectangular discretization of the set 
% [-2.5,2.5]x[-4.5,4.5].This set  is the domain of the Kernel estimator of 
% the invariant density of the SDE defined above. This set it was used be-
% cause the aforementioned estimator converges to the solution of the 
% Fokker-Planck equation.

%% Simulation of the theoretical density: the following function executes 
%  the finite differences method to resolve the Fokker-Planck equation 
[dT,dTx,dTy] = dT4(P,dy,dox,doy);                                                                 
%% Graphics
figure(1)
close all
li=min(min(dT));
ls=max(max(dT));        
zlim([li ls]);
%set(gca,'NextPlot','replaceChildren');
surfl(dox,doy,dT);                                                            
shading interp
colormap gray
title('numerical solution of the Fokker Planck equation')
figure(2)
plot(dox,dTx,'r')
title('numerical solution of the Fokker Planck equation with respect to x')
figure(3)
plot(doy,dTy,'b')
title('numerical solution of the Fokker Planck equation with respect to y')
%% References: 
% Title: Solution of Fokker-Planck equation by finite element and finite 
%        difference methods for nonlinear systems
% Journal: Sadhana
% Year: 2006
% Volume: 31
% This  article can be viewed in mendeley.
toc                                                                        %% returns the elapsed time in seconds of this function (Outputs_FluctuationDissipation.m)
ElapsedTime=toc/60;  
save('DensityResults','dT','dTx','dTy','P','ElapsedTime')                  %%saved variables 
end
