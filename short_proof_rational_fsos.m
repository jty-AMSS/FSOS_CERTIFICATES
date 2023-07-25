function [U,V,S,T,L1normVictor]=short_proof_rational_fsos(f,d,range,k,opti)
if nargin<=3
    k=length(f);
end
if nargin<=4
    opti=0;
end
f=1*f-opti+1/2;
range=range-opti+1/2;
range=range(range>0);
T=select_T(f,d,range,k);
suppf=find(f);
S= SetAdd(suppf,T,f.n);
S=unique(S,'rows');
[U,V,cvx_optval,L1normVictor]=GenerateConidition_rational_complex(f,S,T);
while  cvx_optval~=1
    k=round(k*1.5);
    T=select_T(f,d,range,k);
    
    S= SetAdd(suppf,T,f.n);
    S=unique(S,'rows');
    [U,V,cvx_optval,L1normVictor]=GenerateConidition_rational_complex(f,S,T);
end
end

function T=select_T(f,d,range,k)
p=range(:);
Poly=Polyfit_inf_norm(p,d);
g=Poly(2)*f+Poly(1);
fd=f;
for i=3:length(Poly)
    fd=fd.*f;
    g=Poly(i)*fd+g;
end
[a,b]=find(g);
[~, sortb_Index]=sort(abs(b),'descend');
k=min(k,size(a,1));
T= a(sortb_Index(1:k),:);
T=unique(T,'rows');
end



function Poly=Polyfit_inf_norm(x,d)
A=[x.^[0:d],-ones(length(x),1);-x.^[0:d],-ones(length(x),1)];
Poly=linprog([zeros(1,d+1),1],A,[x.^0.5;-x.^0.5]);
Poly=Poly(1:d+1);
if isnan(Poly(1))&&d==1
    Poly(1)=0;
    Poly(2)=1;
end
end