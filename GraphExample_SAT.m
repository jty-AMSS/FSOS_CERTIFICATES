
Table=[];
Tabler=[];
for n=9:2:49
    G=zeros(n);
    for i=1:n
        for j=1:n
            if (i~=j)&&(abs(mod(i-j,n))<2)
                G(i,j)=1;
            end
        end
    end
    G(n+1,n+1)=0;
    G(n+1,1:n)=1;
    G=G+G';
    G(G>0)=1;
    G=tril(G);
    [Ei,Ej]=find(G);
    Index=1:3*(n+1);
    Index=reshape(Index,[n+1,3]);
    D={};
    for k=1:length(Ei)
        for j=1:3
            D{end+1}=[-Index(Ei(k),j),-Index(Ej(k),j)];
        end
    end
    for k=1:n+1
        D{end+1}=[Index(k,1),Index(k,2),Index(k,3)];
    end
    for j=1:(n+1)
        for k=2:3
            D{end+1}=[-Index(j,k-1),-Index(j,k)];
        end
    end
        
    File=['temp_color.wcnf'];
    fid = fopen(File,'w+');
    fprintf(fid,'c Weighted CNF');
    fprintf(fid,'\n');
    fprintf(fid,'c from my generator');
    fprintf(fid,'\n');
    fprintf(fid,'p wcnf  ');
    fprintf(fid,num2str(3*n+3));
    fprintf(fid,'  ');
    fprintf(fid,num2str(length(D)));
    for i=1:length(D)
        fprintf(fid,'\n');
        y=[1,D{i},0];
        fprintf(fid,num2str(y));
    end
    fclose(fid);
    
    
    

    
    w=ones(length(D),1);
 
    
    f=CNFCharacterFunction2(D,3*(n+1),w);
    S=0:sum(w);
    tic;
    [Q,lb,Index,A,B,g,Suppf_lb,X0]=short_proof_poly_fsos(f,2,S,4,1);
    t=toc;
   Table=[Table;n,t,lb,length(Q)]
    
    
    
    %% 旧的算法
%     [Q,lb,Index,A,B,~,Suppf_lb,X0]=FastLowerBound_Poly(f,[2],S,m,4*m);
%         %     [Q,lb,Index,A,B,g,Suppf_lb,X0]=short_proof_poly_fsos(f,2,S,4);
%         t=toc;
%         Table=[Table;n,t,lb,length(Q)];
%         disp(n);
    %% 新的算法

% [Index,f,w]=Basis08(File,'p');
% tic;
%  [P,A,b]=ComputeSOSByCVX(f,Index);
%  toc;
%  lb=f(f.n*0)-trace(P);
% Table=[Table;n,t,lb,length(P)];
end