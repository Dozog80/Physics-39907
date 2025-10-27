function [L,D,U] = ldu(A)
% ldu <LDU decomposition of a matrix.>
% Usage:: [L,D,U]=ldu(A)
%
% revision history:
% 10/26/2023 Mark D. Shattuck <mds> ldu.m

%% Main
[L,U] = lu(A); % Returns the L and U so A = LU
D = diag(diag(U)); % Constructing D from U, D is a diagonal matrix with 
%entries from the diagonal entries of U
U = D \ U; % Redefine U now so that UD = the old U so LDU = LU_old = A
end
