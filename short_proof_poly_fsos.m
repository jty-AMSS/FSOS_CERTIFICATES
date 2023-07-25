function [Q,lb,Index,A,B,g,Suppf_lb,X0]=short_proof_poly_fsos(f,d,range,ADMmaxiterCoeff,opti,kstart)
if ischar(f)
    [f,w,~]=DICMS2function(f);
    range=0:sum(w);
end
if nargin>=5
    if opti~=0
        f=1*f-opti;
        range=range-opti;
    end
end
f=f+0.5;
range=range(range>=0);
range=range(:);
range=range+0.5;
if nargin==6
    k=kstart;
    
else
    k=length(f);
end
if ADMmaxiterCoeff<10
    ADMmaxiter=ADMmaxiterCoeff*k;
else
    ADMmaxiter=ADMmaxiterCoeff;
end
[Q,lb,Index,A,B,g,Suppf_lb,X0]=FastLowerBound_Poly(f,d,range,k,ADMmaxiter);
while lb<-0.5
    k=round(k*1.5);
    ADMmaxiter=ADMmaxiterCoeff*k;
    [Q,lb,Index,A,B,g,Suppf_lb,X0]=FastLowerBound_Poly(f,d,range,k,ADMmaxiter);
end
end