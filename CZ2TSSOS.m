
function [str,f0]=CZ2TSSOS(f)
% transform C[G] to function which can be computed by TSSOS
N=f.n;

n=length(N);

%% Z_2^n
if all(N==2)
    f=sym(f);
    str=char(vpa(f));
    x=sym('x',[1,n]);
    for i=n:-1:1
        str= strrep(str,char(x(i)), ['x[' num2str(i),']']);
    end
    str= strrep(str,'i','*im');
    return
end
f=sym(f);
f=f+conj(f);
f=1/2*f;





x=sym('x',[1,2*n]);
f=subs(f,conj(x(1:n)),x(n+1:2*n));




ff=sym(f);
x=sym('x',[2,2*n,2],'real');
x=x(1,:,1);
varx=sym('x',[1,2*n]);
f0=subs(ff,varx(n+1:2*n),x(1:n))+subs(ff,varx(n+1:2*n),x(n+1:2*n));
f0=f0/2;
disp('f=')
strf=char(vpa(f0,10));
strf=strrep(strf,'x1_','x[');
strf=strrep(strf,'_1',']');
disp(strf)



str=char(vpa(f,10));
for i=2*n:-1:1
    str= strrep(str,char(x(i)), ['x[' num2str(i),']']);
end
str= strrep(str,'i','*im');

for i=1:n
    fprintf('g%i=2-x[%i]^%i- x[%i]^%i\n',[i,i,N(i),n+i,N(i)]);
    fprintf('h%i=im*x[%i]^%i- im*x[%i]^%i\n',[i,i,N(i),n+i,N(i)]);
end
fprintf('pop=[f')
for i=1:n

fprintf(',g%i,h%i',[i,i])

end
fprintf(']')
fprintf('\n')
    
    

end