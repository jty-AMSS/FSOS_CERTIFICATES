function [g,A,b]=Grid(Q,A,b,n)
%min tr(Q)+||AQ(:)-b)||_1-len(Q)*min_eig(Q)
% [a,~]=eig(Q);
if nargin==3
n=length(Q);
end
y=A*Q(:)-b;
 y2=A'*soft_sign(y,Q);
%% funny random:
% if randi(3)>2
%     y2=A'*soft_sign(y,Q);sign(y);
% else
%    if randn()>0
%        y2=A'*sign(y);
%    else
%        y2=A'*hard_sign(y,Q);
%    end
% end
%% funny end
[a,mineig]= eig(Q,'vector');
[~,k]=min(mineig);
tol=5e1;
%% fun start
tol=max(mineig)-min(mineig);
tol=tol/2;
tol=max(tol,1);
%% fun end

u=a(:,abs(mineig-mineig(k))<tol);
mineig=mineig(abs(mineig-mineig(k))<tol);
mineig=1-abs(mineig-mineig(k));
%% changed
% u=u*abs(mineig(:)+eps);%
% u=u*diag(abs(mineig(:)+eps).^10);

[~,mineig]= eig(Q,'vector');
mineig=mineig(abs(mineig-mineig(k))<tol);
mineig=mineig-mineig(k);
mineig=-1*mineig;
u=u*diag(0.2.^abs(mineig(:)+eps));
%% changed ended
u=u*u';
u=u/(norm(u(:)));
G1=-n*u;
n=length(Q);
G2=reshape(y2,[n,n]);
G3=eye(n);
% G3(1)=0;
g=G1+G2+G3;
g=-g;
% g=-g;

end

function z=soft_sign(y,Q)
tol=1/length(Q)^2;
y1=y;
y(abs(y)>tol)=0;
y1=y1-y;
y1=y1/tol;
z=sign(y1)+y;
end

function z=hard_sign(y,Q)
tol=1/length(Q)^2;
y(abs(y)<tol)=0;
z=sign(y);
end


