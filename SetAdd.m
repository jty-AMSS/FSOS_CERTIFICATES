
function S=SetAdd(A,B,N)
%returns S={a+b;(a,b)\in A*B} ,A,B \subset prod(Z_N(i))
n=length(N);
na=size(A,1);
nb=size(B,1);
S=kron(A,ones(nb,1))+kron(ones(na,1),B);
S=mod(S,N);
end