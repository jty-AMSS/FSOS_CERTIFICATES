function [Index,A,B]=FastLowerBound_Poly_Index(f,d,range,k)
%Êä³öÎªIndex£¬A,B
%range: range of function £¬class£º 1-d array
% out put the first k Index of P(f)
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
Index= a(sortb_Index(1:k),:);
[A,B,~,~,Suppf_lb]=GenerateConidition(f,Index);
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


function [A,b,Eq,Suppf_Minus_Condition,Suppf_lb]=GenerateConidition(f,Index,ifconsider_complex)
    Eq=containers.Map('KeyType',  'char', 'ValueType', 'any');
    m=size(Index,1);
    N=f.n;
    for i=1:m
        for j=1:m
            t=mod(Index(i,:)-Index(j,:),N);
            if isKey(Eq,char(t))
                Eq(char(t))=[Eq(char(t)),m*(i-1)+j];
            else
                Eq(char(t))=[m*(i-1)+j];
            end
        end
    end
    suppf=find(f);
    G=keys(Eq);
%     %%%
    Eq.remove(char(zeros(1,length(f.n))));
%     %%%
    G=keys(Eq);
    Ax=[];Ay=[];b=zeros(length(G),1);
    for i=1:length(G)
        t=G{i};
        Ax(end+1:end+(length(Eq(t))))=i;
        Ay(end+1:end+(length(Eq(t))))=Eq(t);
        b(i)=f(t);
    end
    A=sparse(Ax,Ay,1,length(G),m*m);
    G{end+1}=char(zeros(1,length(f.n)));
    Suppf_Minus_Condition=setdiff(suppf,cell2mat(G(:))+0,'rows');
    Suppf_lb=0;
    for i=1:size(Suppf_Minus_Condition,1)
        Suppf_lb=Suppf_lb+abs(f(char(Suppf_Minus_Condition(i,:))));
    end
    return
end

