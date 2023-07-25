diary on
for n=10:10:100
    N=ones(1,n)+1;
    f=CZ(N);
    g=CZ(N);
    for i=1:n
        e=g.n;e(i)=1;
        g(e)=1/2;
    end
    f=(1*g)*(1*g-1);
    f=f+1/2;
    range=0:(n/2*(n/2+1));
    range=range+1/2;
    [Q,lb(n/10),Index,A,B]=FastLowerBound_Poly(f,2,range,n*n/2,1e4);
    disp(lb');
end
for n=10:10:100
    N=ones(1,n)+1;
    f=CZ(N);
    g=CZ(N);
    for i=1:n
        e=g.n;e(i)=1;
        g(e)=1/2;
    end
    f=(1*g)*(1*g-1);
    range=0:(n^2/4 + n/2);
    f=f+1/2;
    range=range+1/2;
    [Q,lb3(n/10),Index,A,B]=FastLowerBound_Poly(f,3,range,n*5,1e4);
    disp(lb3');
end
diary off
