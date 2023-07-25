function [f,S]=Coloring2GpFun(G,k)
%最大子图k染色问题转化为 群函数f，具体是：
%\delta (x,y): C_k*C_k->Z , \delta (x,y)=1 <=> x==y and \delta (x,y)=0 <=> x!=y
% f=\sum(u,v \in E(G)) \delta(u,v);

%input:  
% G: (sparse) matrix, which is the adjacency matrix
% k: integer, coloring number

% output:
% f: C[C_k^{|V|}]->Z, with f(g)= the number of deleting with coloring
% method g,
% S: im(f) \subset S
G= tril(G);
[edge_i,edge_j]=find(G);
N=ones(1,length(G));
N=k*N;
S=(0:sum(G(:)));
n=length(G);
f=CZ(N);
for  i=1:length(edge_i)
    u=edge_i(i);v=edge_j(i);
for j=1:k
    t=zeros(1,n);
    t(u)=j;
    t(v)=k-j;
    f(t)=f(t)+1/k;
end
end
end