
Table=[];
Tabler=[];
for n=5:25
    G=zeros(n);
    G=G+1;
    G=G-diag(diag(G));
    [f,S]=Coloring2GpFun(G,n-1);
    f=f-1;
    f=f+0.5;
    S=S-0.5;
    S=S(S>0);
    tic;
    m=round(3*length(f));
    [Q,lb,Index,A,B,g,Suppf_lb,X0]=short_proof_poly_fsos(f,2,S,4);
    t=toc;
    Table=[Table;n,t,lb,length(Q)];
end