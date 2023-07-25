function [err,p,q,c,g,HH,G]=CheckFSOS_rational(f,S,T,U,V)

Q=U;Index=S;
n=length(Q);
if min(eig(Q))<0
    Q=Q-min(eig(Q))*eye(n);
end
H=chol(Q);

g=cell(1,n);
fsos=CZ(f.n);
for i=1:n
    tempg=CZ(f.n);
    tempg2=CZ(f.n);
    for k=1:size(Index,1)
        tempg(char(Index(k,:)))=H(i,k);
    end
    g{i}=tempg;
    tempg2=tempg*conj(tempg);
    fsos=fsos+tempg2;
end
p=fsos;
HH=g;
g={};
Q=V;Index=T;
n=length(Q);
if min(eig(Q))<0
    Q=Q-min(eig(Q))*eye(n);
end
H=chol(Q);
g=cell(1,n);
fsos=CZ(f.n);
for i=1:n
    tempg=CZ(f.n);
    tempg2=CZ(f.n);
    for k=1:size(Index,1)
        tempg(char(Index(k,:)))=H(i,k);
    end
    g{i}=tempg;
    tempg2=tempg*conj(tempg);
    fsos=fsos+tempg2;
end
G=g;
q=fsos;
err=p-q*f;

c=coeffs(sym(err));
end