function [P,Index,A,B]=short_proof_poly(f,d,range,k,L)
if nargin==4
    L=0;
end
f=1*f;
f=f-(L-1/2);
range=range-(L-1/2);
range=range(range>0);
[Index,A,B]=FastLowerBound_Poly_Index(f,d,range,k);
M=size(Index,1);
% cvx_begin
% variable  P(M,M) hermitian semidefinite
% minimize(1)
% norm(A*P(:)-B,1)<=1/2;
% cvx_end

%% Opt Prob:
t=length(f.n);
t=zeros(1,t);
f0=f(t);
cvx_begin
variable  P(M,M) hermitian semidefinite
maximize(f0-trace(P)-norm(A*P(:)-B,1))
cvx_end

end