% Block upper triangularization of 3 matrices
% Devoted to the Draft Ha18_difference
% function [r2,r1,r0,A2,A1,A0] = BlkUpTri(A2,A1,A0) bring the systems to
% the block upper triangular form
% NOTICE: 3 matrices A2,A1,A0 must of of the same size
% BlkUptri(A2,A1,A0) returns the block upper triangular form
% [A21   A11   A01]    r2 rows
% [0     A12   A02]    r1 rows
% [0     0     A03]    r0 rows
% [0     0     0  ]    v rows
% The matrices on the main diagonals have full row rank.
% Copyright: @Phi Ha

function [r2,r1,r0,A2,A1,A0] = BlkUpTri(A2,A1,A0)

r2 = rank(A2); 
Z11 = orth(A2); Z12 = null(A2'); Z1 = [Z11 Z12];
A2 = Z1' * A2; 
A1 = Z1' * A1; 
A0 = Z1' * A0;

A12 = A1(r2+1:end,:); r1 = rank(A12);
Z21 = orth(A12); Z22 = null(A12'); Z2 = [Z21 Z22]; 
A2(r2+1:end,:) = Z2' * A2(r2+1:end,:);
A1(r2+1:end,:) = Z2' * A1(r2+1:end,:);
A0(r2+1:end,:) = Z2' * A0(r2+1:end,:);

A03 = A0(r2+r1+1:end,:); r0 = rank(A03);
Z31 = orth(A03); Z32 = null(A03'); Z3 = [Z31 Z32];
A2(r2+r1+1:end,:) = Z3' * A2(r2+r1+1:end,:);
A1(r2+r1+1:end,:) = Z3' * A1(r2+r1+1:end,:);
A0(r2+r1+1:end,:) = Z3' * A0(r2+r1+1:end,:);


