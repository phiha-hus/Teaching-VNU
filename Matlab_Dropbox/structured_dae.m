clear all; close all; clc

syms c m k g

E = [c m 0; m 0 0;0 0 0]

A = [-k 0 g; 0 -m 0; g 0 0]

E*A

A*E

simplify(E*A-A*E)