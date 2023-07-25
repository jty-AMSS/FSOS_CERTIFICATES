function [U,V,cvx_optval,L1normVictor]=GenerateConidition_rational_complex(f,S,T)
%output:
%A(:)*U(:)=B(:)*V(:)
N=f.n;
suppf=find(f);
SmS=SetAdd(S,-S,N);
TmT=SetAdd(T,-T,N);
TmTpSuppf=SetAdd(TmT,suppf,N); 

lambda_Index=[SmS;TmTpSuppf];
lambda_Index=unique(lambda_Index,'rows');

Equi=containers.Map('KeyType',  'char', 'ValueType', 'any');
Equj=containers.Map('KeyType',  'char', 'ValueType', 'any');
Eqvi=containers.Map('KeyType',  'char', 'ValueType', 'any');
Eqvj=containers.Map('KeyType',  'char', 'ValueType', 'any');
Eqvconst=containers.Map('KeyType',  'char', 'ValueType', 'any');
for i=1:size(S,1)
    for j=1:size(S,1)
        alpha=S(i,:);
        beta=S(j,:);
        t=mod(alpha-beta,N);
        if isKey(Equi,char(t))
             Equi(char(t))=[Equi(char(t)),i];
             Equj(char(t))=[Equj(char(t)),j];
        else
            Equj(char(t))=[j];  Equi(char(t))=[i]; 
        end
    end
end

for i=1:size(T,1)
    for j=1:size(T,1)
        for k=1:size(suppf,1)
        alpha=T(i,:);
        beta=T(j,:);
        gamma=suppf(k,:);
        t=mod(alpha-beta+gamma,N);
        if isKey(Eqvi,char(t))
             Eqvi(char(t))=[Eqvi(char(t)),i];
             Eqvj(char(t))=[Eqvj(char(t)),j];
             Eqvconst(char(t))=[Eqvconst(char(t)),f(gamma)];
        else
            Eqvj(char(t))=[j];  Eqvi(char(t))=[i];Eqvconst(char(t))=[f(gamma)];
        end
        end
    end
end

cvx_begin
variable  U(size(S,1),size(S,1)) hermitian semidefinite
variable  VV(size(T,1),size(T,1)) hermitian semidefinite
V=VV+eye(size(T,1));
minimize(1)
L1normVictor=zeros(size(lambda_Index,1))*U(1,1);
for i=1:size(lambda_Index,1)
    lambda=lambda_Index(i,:);
    lambda=char(lambda);
    %LHS

    if isKey(Equi,lambda)
        ui=Equi(lambda);
        uj=Equj(lambda);
        lhs=trace(U(ui,uj));
    else
        lhs=0;
    end
    %RHS
    if isKey(Eqvi,lambda)
        vi=Eqvi(lambda);
        vj=Eqvj(lambda);
        vc=Eqvconst(lambda);
        rhs=trace(diag(vc)*V(vi,vj));
    else
        rhs=0;
    end
    L1normVictor(i)=(lhs-rhs);
    
end
norm(L1normVictor,1)<=0.5;
cvx_end
L1normVictor=norm(L1normVictor,1)
% [u]=CheckFSOS(CZ(f.n),S,U);u=-1*u;
% [v]=CheckFSOS(CZ(f.n),T,V);v=-1*v;
% err= v*f-u;errc=coeffs(sym(err));vpa(norm(errc,1))

% % m=size(Index,1);
% % N=f.n;
% % for i=1:m
% %     for j=1:m
% %         t=mod(Index(i,:)-Index(j,:),N);
% %         if isKey(Eqi,char(t))
% %             Eqj(char(t))=[Eqj(char(t)),i];
% %             Eqi(char(t))=[Eqi(char(t)),j];
% %         else
% %             Eqj(char(t))=[i];  Eqi(char(t))=[j]; %(j,i)->m*(i-1)+j
% %         end
% %     end
% % end
% % 
% % suppf=find(f);
% % G=keys(Eqi);
% % %%%
% % Eqi.remove(char(zeros(1,length(f.n))));
% % 
% % Eqj.remove(char(zeros(1,length(f.n))));
% % %%%
% % G=keys(Eqi);
% % Ax=[];Ay=[];b=zeros(2*length(G),1);
% % Axm=[];Aym=[];
% % G_m=length(G);
% % for i=1:length(G)
% %     t=G{i};
% %     Ax(end+1:end+(length(Eqi(t))))=i;
% %     Ay(end+1:end+(length(Eqi(t))))=Eqi(t)+2*m*(Eqj(t)-1); %(i,j)_2m
% %     Ax(end+1:end+(length(Eqi(t))))=i;
% %     Ay(end+1:end+(length(Eqi(t))))=m+Eqi(t)+2*m*(Eqj(t)-1+m);%(i+m,j+m)
% %     
% %     Ax(end+1:end+(length(Eqi(t))))=i;
% %     Ay(end+1:end+(length(Eqi(t))))=m+Eqj(t)+2*m*(Eqi(t)-1+m);%(j+m,i+m)
% %     Ax(end+1:end+(length(Eqi(t))))=i;
% %     Ay(end+1:end+(length(Eqi(t))))=Eqj(t)+2*m*(Eqi(t)-1);%(j,i)
% %     
% % % % % % % % % % % % %     
% %     Axm(end+1:end+(length(Eqi(t))))=length(G)+i;
% %     Aym(end+1:end+(length(Eqi(t))))=m+Eqi(t)+2*m*(Eqj(t)-1);%(i+m,j)
% %     Ax(end+1:end+(length(Eqi(t))))=length(G)+i;
% %     Ay(end+1:end+(length(Eqi(t))))=Eqi(t)+2*m*(Eqj(t)-1+m);%(i,m+j)
% %     
% %     Axm(end+1:end+(length(Eqi(t))))=length(G)+i;
% %     Aym(end+1:end+(length(Eqi(t))))=Eqj(t)+2*m*(m+Eqi(t)-1);%(j,i+m)
% %     Ax(end+1:end+(length(Eqi(t))))=length(G)+i;
% %     Ay(end+1:end+(length(Eqi(t))))=m+Eqj(t)+2*m*(Eqi(t)-1);%(m+j,i)
% %     b(i)=real(f(t));
% %     b(G_m+i)=imag(f(t));
% % end
% % Ap=sparse(Ax,Ay,1/2,2*length(G),4*m*m);
% % Am=sparse(Axm,Aym,-1/2,2*length(G),4*m*m);
% % A=Ap+Am;
% % G{end+1}=char(zeros(1,length(f.n)));
% % Suppf_Minus_Condition=setdiff(suppf,cell2mat(G(:))+0,'rows');
% % Suppf_lb=0;
% % for i=1:size(Suppf_Minus_Condition,1)
% %     Suppf_lb=Suppf_lb+abs(f(char(Suppf_Minus_Condition(i,:))));
% % end
end


function S=SetAdd(A,B,N)
%returns S={a+b;(a,b)\in A*B} ,A,B \subset prod(Z_N(i))
n=length(N);
na=size(A,1);
nb=size(B,1);
S=kron(A,ones(nb,1))+kron(ones(na,1),B);
S=mod(S,N);
end