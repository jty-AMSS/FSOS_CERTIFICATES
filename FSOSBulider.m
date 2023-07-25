function [lb_l,Q,f,X0,Index,x,y]=FSOSBulider(path,d,sparsity,opti,ADMmaxiter)
if nargin==3
    opti=0;
end

[f,w,~]=DICMS2function(path);
if (sparsity<3)&&(sparsity>0)
    sparsity=sparsity-1e-10;
end
if abs(sparsity-round(sparsity))>0
    sparsity=length(f)*sparsity;
    sparsity=round(sparsity);
end
if sparsity==-1
    sparsity=length(f);
end
if nargin<5
    ADMmaxiter=1.5*sparsity;
end

if opti~=0
f=1*f-opti;
end
f=f+0.5;
range=0:sum(w);
range=range-opti;
range=range(range>=0);
range=range(:);
range=range+0.5;
[Q,~,Index,A,B,~,Suppf_lb,X0]=FastLowerBound_Poly(f,d,range,sparsity,ADMmaxiter);
lb_l=norm(A*Q(:)-B,1)+trace(Q)- length(Q)*lambda_min(Q);
f_const=f((zeros(1,length(f.n))));
lb_l=f_const-lb_l;
lb_l=lb_l-Suppf_lb;
end