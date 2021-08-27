%% Introduction
% David Vadnais
% extended Homework 1
% Problem 3
clear; clc;

%% Part A
% GIVENS
% Strain CWID = ****2326
c11 = .4450;
c12 = .3950;
c13=.4050;
c33=.4440;
c44=0.0655;
c66=0.1220;
%tetragonal solid
c22= c11;
c23=c13;
c55=c44;

ciijj = (c11+c12+c13+c22+c23+c13+c23+c33)*10^10
cijij= (c11+c66+c55+c66+c22+c44+c55+c44+c33)*10^10

lambV = 1/15*(2*ciijj-cijij)
muV=1/30*(3*cijij-ciijj)

disp('Vogt average elastic modulus (GPa)')
disp(lambV/10^9)
disp('Vogt average rigidity modulus (GPa)')
disp(muV/10^9)

%% Part B
c = [c11,c12,c13,0,0,0;...
    c12,c22,c23,0,0,0;...
    c13,c23,c33,0,0,0;...
    0,0,0,c44,0,0;...
    0,0,0,0,c55,0;...
    0,0,0,0,0,c66]*10^10;
s = inv(c)

siijj = s(1,1)+s(1,2)+s(1,3)+s(1,2)+s(2,2)+s(2,3)+s(1,3)+s(2,3)+s(3,3)
sijij = s(1,1)+s(6,6)+s(5,5)+s(6,6)+s(2,2)+s(4,4)+s(5,5)+s(4,4)+s(3,3)

EReuss = 1/((1/15)*(2*sijij+siijj))
GReuss = 1/((1/15)*(6*sijij-2*siijj))
disp('Reuss average elastic modulus (GPa)')
disp(EReuss/10^9)
disp('Reuss average rigidity modulus (GPa)')
disp(GReuss/10^9)

%% Part C
eps = [823,1023,-1023;1023,526,-923;-1023,-923,-526]*10^-6;%strains

epsV = [eps(1,1);eps(2,2);eps(3,3);2*eps(2,3);2*eps(1,3);2*eps(1,2)];

Sigma = c*epsV;

TensSigma = [Sigma(1),Sigma(6),Sigma(5);...
            Sigma(6),Sigma(2),Sigma(4);...
            Sigma(5),Sigma(4),Sigma(3)];
%Invariants
I1 = trace(TensSigma);
I2 = TensSigma(1,1)*TensSigma(2,2)+TensSigma(1,1)*TensSigma(3,3)+...
    TensSigma(2,2)*TensSigma(3,3)-TensSigma(1,2)^2-TensSigma(1,3)^2-...
    TensSigma(2,3)^2;
I3 = det(TensSigma);

sigP = roots([1,-I1,I2,-I3]);

[Vec,Val]=eig(TensSigma)

sig1 = Val(3,3);
sig2 = Val(2,2);
sig3 = Val(1,1);

disp('prinicpal stress 1 (MPa)')
disp(sig1/10^6)
disp('prinicpal stress 2 (MPa)')
disp(sig2/10^6)
disp('prinicpal stress 3(MPa)')
disp(sig3/10^6)



%% Part D
v=(EReuss-GReuss)/(2*GReuss)

S=[1/EReuss,-v/EReuss,-v/EReuss,0,0,0;...
       -v/EReuss,   1/EReuss,-v/EReuss,0,0,0;...
       -v/EReuss,-v/EReuss,   1/EReuss,0,0,0;...
       0,0,0,1/GReuss,0,0;...
       0,0,0,0,1/GReuss,0;...
       0,0,0,0,0,1/GReuss]
   
C = inv(S);

SigmaD = C*epsV

TensSigma = [SigmaD(1),SigmaD(6),SigmaD(5);...
            SigmaD(6),SigmaD(2),SigmaD(4);...
            SigmaD(5),SigmaD(4),SigmaD(3)];

% stress invariants
I1D = trace(TensSigma);
I2D = TensSigma(1,1)*TensSigma(2,2)+TensSigma(1,1)*TensSigma(3,3)+...
    TensSigma(2,2)*TensSigma(3,3)-TensSigma(1,2)^2-TensSigma(1,3)^2-...
    TensSigma(2,3)^2;
I3D = det(TensSigma);
sigP = roots([1,-I1,I2,-I3])

disp('prinicpal stress 1 (MPa)')
disp(sigP(1)/10^6)
disp('prinicpal stress 2 (MPa)')
disp(sigP(2)/10^6)
disp('prinicpal stress 3(MPa)')
disp(sigP(3)/10^6)




