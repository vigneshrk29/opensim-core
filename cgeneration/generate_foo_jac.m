% This script generates c-code from the function foo.m that contains the 
% expression graph of function F. The c-code foo_jac will contain the 
% function F and its Jacobian in a format that can be exploited by CasADi.
% Authors: Joris Gillis and Antoine Falisse

function [] = generate_foo_jac(dim)

import casadi.*
cg = CodeGenerator('foo_jac');
% arg should have the dimensions of the combined inputs of F, i.e. NX + NU
arg = SX.sym('arg',dim); 
[y,~,~] = foo(arg);
F = Function('F',{arg},{y});
cg.add(F);
cg.add(F.jacobian())
cg.generate();
