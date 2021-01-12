% This script generates c-code from the function foo.m that contains the 
% expression graph of function F. The c-code foo_jac will contain the 
% function F and its Jacobian in a format that can be exploited by CasADi.
% Authors: Joris Gillis and Antoine Falisse

% function [] = generate_foo_jac(dim)

dim = 18;

import casadi.*
cg = CodeGenerator('foo_jac');
% arg should have the dimensions of the combined inputs of F, i.e. NX + NU
arg = SX.sym('arg',dim); 
[y_T,~,b_T] = foo_T(arg);
[y_F,~,b_F] = foo_F(arg);
% F = Function('F',{arg},{y});
% cg.add(F);
% cg.add(F.jacobian())
% cg.generate();

% y = if_else(b_T==1, y_T, y_F);
F_T = Function('F',{arg},{y_T});
F_F = Function('F',{arg},{y_F});

% cg.add(F);
% cg.add(F.jacobian())
% cg.generate();

angle = 250;
in1 = [cosd(angle), 0, sind(angle), 0, 1, 0, -sind(angle), 0, cosd(angle)];
in2 = [cosd(-angle), 0, sind(-angle), 0, 1, 0, -sind(-angle), 0, cosd(-angle)];
out_T = full(F_T([in1, in2]));


angle = 20;
in3 = [cosd(angle), 0, sind(angle), 0, 1, 0, -sind(angle), 0, cosd(angle)];
in4 = [cosd(-angle), 0, sind(-angle), 0, 1, 0, -sind(-angle), 0, cosd(-angle)];
out_F = full(F_F([in3, in4]));  
out_F_250 = full(F_F([in1, in2]));  

% out_T_250 = full(F_T([in1, in2]));  
% out_F_250 = full(F_F([in1, in2]));
% out_250_in_20 = full(F_250_in([in3, in4]));
% out_250_out_20 = full(F_250_out([in3, in4]));


% 
% assert(round(out1, 5)==round(2.44346,5))
% assert(round(out2, 5)==round(0.698132,5))
% 

