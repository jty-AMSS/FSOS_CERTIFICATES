function  c=TSSOS_Code(f,d)
% ‰»Î£∫ f \in CZ[N]
%  ‰≥ˆ£∫ TSSOS¥˙¬Î
N=f.n;
n=length(N);
N=N(:);
n=n*2;
x=sym('x',[2,n,2],'real');
x=x(1,:,1);
y=reshape(x,[n/2,2]);
y=y*[1;1i];
C=expand(y.^N-1);
C=[real(C);imag(C)];
strC=char(C(:).');
strC=strC(9:end-2);
strC=strrep(strrep(strC,'x1_','x['),'_1',']');
f=sym(f);
strx=sym('x',[n/2,1],'real');
f=subs(f,strx,y);
f=expand(f+conj(f));
f=f/2;
strf=strrep(strrep(char(vpa(f,10)),'x1_','x['),'_1',']');
c={};
c{1}='using TSSOS';
c{2}='using DynamicPolynomials';
c{3}=['n=' num2str(n)];
c{4}='@polyvar x[1:n]';
c{5}=['f=' strf];
c{6}=['C=' strC];
c{7}=['pop=[f;C]'];
c{8}=['neq=' num2str(length(C))];
c{9}=['d=' num2str(d)];
c{10}='import Dates;x1=time();';
c{11}='opt,sol,data = tssos_first(pop,x, d, numeq=neq,TS="MD",Gram=true)';
c{12}='x2=time();';
c{13}='T21=x2-x1;;basis=data.basis;sparsity=size(basis[1]);sparsity=convert(Float64,sparsity[2]);print([T21,opt,sparsity])';
c{14}=[' '];
for i=1:length(c)
    disp(c{i});
end
end