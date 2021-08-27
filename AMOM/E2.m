%% Introduction
% David Vadnais
% extended Homework 1
clear; clc;

%% Part A
% GIVENS
% CWID = ****2326
A=23; B=26;
sigma = [-A,-15,200+A;-15,B,A;200+A,A,20];%Stress Tensor
n = transpose([1/sqrt(3),1/sqrt(3),1/sqrt(3)]);

% TRACTION
t = sigma*n%MPa

%% Part B - Normal stress
sigmaN = norm(t.*n)

%% PArt C - Shear Stress
tau = sqrt(norm(t)^2-sigmaN^2)

%Example problem from notes
% sigma = [-15,0,-25;0,10,-20;-25,-20,0]
% n=[0,4/5,3/5]'
% t= sigma*n
% sigmaN=sum(t.*n)
% tau = sqrt(norm(t)^2-sigmaN^2)