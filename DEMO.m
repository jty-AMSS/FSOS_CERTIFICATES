
%% Example 3.5
N=[2,2,2];
f=CZ(N);
f([0,0,0])=13/8;
f([1,0,0])=3/8;
f([0,1,0])=3/8;
f([0,0,1])=3/8;
f([1,0,1])=1/8;
f([1,1,0])=1/8;
f([0,1,1])=1/8;
f([1,1,1])=-1/8;
S=1:3;
[Q,lb,Index,A,B]=short_proof_poly_fsos(f,1,S,4,1,4);
[errf,F]=CheckFSOS(f-1,Index,Q);

%% Example 3.15
N=ones(1,6)+1;
g=CZ(N);
E=eye(6);
g(E(1,:))=1/4;
g(E(2,:))=1/4;
for k=3:6
    g(E(k,:))=1/2;
end
g(E(1,:)+E(2,:))=1/4;
g(g.n*0)=3.25;
[Q,lb,Index,A,B]=short_proof_poly_fsos(g,1,S,4,1,5);
[errg,G]=CheckFSOS(1*g-1,Index,Q);

%% Example 3.5
N=[2,2,2];
f=CZ(N);
f([0,0,0])=13/8;
f([1,0,0])=3/8;
f([0,1,0])=3/8;
f([0,0,1])=3/8;
f([1,0,1])=1/8;
f([1,1,0])=1/8;
f([0,1,1])=1/8;
f([1,1,1])=-1/8;
% S=1:3;
% [U,V,S,T,L1normVictor]=short_proof_rational_fsos(f,2,S,4,1);
% [errpq,p,q,c]=CheckFSOS_rational(f,S,T,U,V);


%% Newman rational


f0=1*f-1/2;
V=CZifft(f0);
m=max(V(:))
m=sym(m);

d=1;
d=d*2;
xi=exp(-(sym(d))^(-1/2));
syms t;
pd=1;
for k=0:(d-1)
    pd=pd*(t+xi^k);
end

pd1=t*(pd-subs(pd,t,-t));pd1=expand(pd1);
pd2=(pd+subs(pd,t,-t));pd2=expand(pd2);

pd1=subs(pd1,t,t^(1/2));
pd2=subs(pd2,t,t^(1/2));

r=sqrt(m)*subs(pd1,t,t/m)/subs(pd2,t,t/m);



p=f0*sqrt(10)*(1 + exp(-sqrt(2)/2));
q=(2*f0 + 5*exp(-sqrt(2)/2));

V=CZifft(q*q);
u=min(V(:));
u=1/u;
errl1=vpa(norm(coeffs(sym(u*((1*f-1/2)*q*q)-u*p*p)),1));

%% Example 5.3
N=[2,2,2,2];
f=CZ(N)
E=eye(4);
for k=1:4
    f(E(k,:))=1;
end
f(E(1,:)+E(2,:))=1;
f(E(1,:)+E(3,:))=2;
f(E(1,:)+E(4,:))=-3;
f(E(1,:)*0)=8;
range=0:18;
L=2;
[U,V,S,T,L1normVictor]=short_proof_rational_fsos(f,1,range,2,L);
[err,p,q,c,g,H,G]=CheckFSOS_rational((1*f-3/2),S,T,U,V);

disp('g_:')
for i=1:length(H)
    disp(vpa(H{i},2))
end
disp('h_:')
for i=1:length(G)
    disp(vpa(G{i},2))
end
disp('error')
disp(vpa(norm(c,1)))