


%% Wheel graph

Table=[];
Tabler=[];
for n=9:2:49
G=zeros(n);
for i=1:n
for j=1:n
if (i~=j)&&(abs(mod(i-j,n))<2)
G(i,j)=1;
end
end
end
G(n+1,n+1)=0;
G(n+1,1:n)=1;
G=G+G';
G(G>0)=1;
[f,S]=Coloring2GpFun(G,3);
f=f-1;
f=f+0.5;
S=S-0.5;
S=S(S>0);
tic;
m=round(3*length(f));
 [Q,lb,Index,A,B,g,Suppf_lb,X0]=short_proof_poly_fsos(f,2,S,1e5);
t=toc;
Table=[Table;n,t,lb,length(Q)];
end