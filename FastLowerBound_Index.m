function [X0,lb]=FastLowerBound_Index(f,Index,ADMmaxiter)
ifconsider_comple=1;
c=coeffs(sym(f));
if all(real(c)==c)&& all(f.n==2)
    ifconsider_comple=0;
end
[A,b]=GenerateConidition(f,Index,ifconsider_comple);
[X0,lb_nal]=ComputeSOSBySDPNAL(A,b,ADMmaxiter);
Q=X0;
n=length(X0);
F=@(x)trace(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')-n*lambda_min(0.5*reshape(x,[n,n])+0.5*reshape(x,[n,n])')+norm(A*x(:)-b,1);
lb= f(zeros(1,length(f.n)))-F(Q);

end