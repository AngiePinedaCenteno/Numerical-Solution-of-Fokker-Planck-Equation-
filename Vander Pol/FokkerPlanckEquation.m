function [dT] = FokkerPlanckEquation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generates a numerically soution of the Fokker-Planck equa-
%tion corresponding to the Van der Pol oscillator:
%        dX_t=Y_tdt
%        dY_t=sigmadB_t-((c1X_t^2-c2)Y_t+w_0^2X_t)dt
%% Input Variables and problem parameters
% - dx,dy:grid spacings for x and y in the discretized domain of density
% - dox: discretized domain in x direction for the density
% - doy: discretized domain in y direction for the density
% - P: vector that contains the parameters of the stochastic differential 
%   equation  defined  above. Specifically, c1=P(1), c2=P(2), sigma=P(3) y
%   W0=P(4).
%% Assignation of the parameters values
P=[0.2,0.2,1,2];
%% Discretization of the  domain of the density. 
dy=1/10;dx=dy;dox=(-4:dy:4);doy=(-7:dy:7); 
%% Remark: 
% To obtain the numerical solution of the Fokker-Planck equation it was 
% used  the finite differences method. The grid used to compute the numeri-
% cal solution  is a rectangular discretization of the set [-4,4]x[-7,7]. 
% This set  is the domain of the Kernel estimator of the invariant density 
% of the SDE defined above. This set it was used because the aforementioned
% estimator converges to the solution of the Fokker-Planck equation.
%% Simulation of the theoretical density
[dT,dTx,dTy] = dT4(P,dy,dox,doy);
%% Graphics
close all
figure(1)
mesh(dox,doy,dT)                                                       
title('numericaly solution of the Fokker Planck equation')
figure(2)
plot(dox,dTx,'r')
grid
title('numericaly solution of the Fokker Planck equation with respect to x')
figure(3)
plot(doy,dTy,'b')
grid
title('numericaly solution of the Fokker Planck equation with respect to y')
end
