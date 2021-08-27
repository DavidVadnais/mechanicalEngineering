%% Introduction
% David Vadnais
% extended Homework 1
% Problem 3
clear; clc;

%% Part A
% GIVENS
% Strain CWID = ****2326
 syms sxx sxy sxz syy syz szz 
syms nxpx nxpy nxpz nypx nypy nypz nzpx nzpy nzpz
% 
T = [nxpx, nxpy, nxpz;...
    nypx, nypy, nypz;...
    nzpx, nzpy, nzpz]
 s = [sxx, sxy,sxz;sxy,syy,syz;sxz,syz,szz]
% 
sprime= T*s*transpose(T)

% syms thetaxpx thetaxpy thetaxpz thetaypx thetaypy thetaypz thetazpx thetazpy thetazpz
% T2 = [cos(thetaxpx), cos(thetaxpy), cos(thetaxpz); ...
%     cos(thetaypx), cos(thetaypy), cos(thetaypz); ...
%     cos(thetazpx), cos(thetazpy),cos(thetazpz)]
% s2 = [sxx, sxy,sxz;sxy,syy,syz;sxz,syz,szz]
% 
% sprime2= T2*s2*transpose(T2)
% sum(diag(sprime2))
% 
% 
% 
% %% 
% syms thetaX thetaY thetaZ
% T3 = [cos(thetaX), sin(thetaX), cos(thetaxpz);...
%     -sin(thetaX), cos(thetaX), cos(thetaypz); ...
%     cos(thetazpx), cos(thetazpy),cos(thetaZ)]
% s3 = [sxx, sxy,sxz;sxy,syy,syz;sxz,syz,szz]
% 
% sprime3= T3*s3*transpose(T3)
% sum(diag(sprime3))

% 
% syms phi theta psi
% 
% Tz = [cos(theta), sin(theta) , 0; -sin(theta) ,cos(theta),0;0,0,1]
% Ty = [cos(phi), 0,sin(phi);0,1,0;-sin(phi),0,cos(phi)]
% Tx = [1,0,0;0,cos(psi),sin(psi);0,-sin(psi),cos(psi)]
% 
% T=Tz*Ty*Tx
%  
% Sp = T*s*transpose(T)
% sumDiagPrim = sum(diag(Sp))
% sumDiag = sum(diag(s))













