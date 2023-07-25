
function [Index,f,w]=Basis08(File,BasisStr)
% BasisStr='GW' ,'p','ap' 't' 'pt'
if nargin<2
    BasisStr='p';
end
[f,w,D]=DICMS2function(File);
Index=PolySOSWithAlpha(f,D,BasisStr);
end
function subs= PolySOSWithAlpha(f,D,BasisStr)
n=f.n;n=length(n);
subs=[];
subs=[subs;char(zeros(1,n))];
subs=[subs;char(eye(n))];
%GW end
if (strcmp(BasisStr,'p')||strcmp(BasisStr,'pt'))
    MpIndex=[];
    Mp=[];
    for i=1:length(D)
        t=D{i};
        if (length(t)>=2)
            t=abs(t);
            t=sort(t);
            MpIndex=[MpIndex;nchoosek(t,2)];
        end
    end
    MpIndex=unique(MpIndex,'rows');
    Mp=zeros(size(MpIndex,1),n);
    for i=1:size(MpIndex,1)
        Mp(i,MpIndex(i,:))=1;
    end
    Mp=char(Mp);
    subs=[subs;Mp];
end
%M_p End
if (strcmp(BasisStr,'t')||strcmp(BasisStr,'pt'))
    MtIndex=[];
    Mt=[];
    for i=1:length(D)
        t=D{i};
        if (length(t)>=3)
            t=abs(t);
            t=sort(t);
            MtIndex=[MtIndex;nchoosek(t,3)];
        end
    end
    MtIndex=unique(MtIndex,'rows');
    Mt=zeros(size(MtIndex,1),n);
    for i=1:size(MtIndex,1)
        Mt(i,MtIndex(i,:))=1;
    end
    Mt=char(Mt);
    subs=[subs;Mt];
end
%M_t End
if (strcmp(BasisStr,'ap'))
    MapIndex=nchoosek(1:(length(f.n)),2);
    for i=1:size(MapIndex,1)
        M_ap(i,MapIndex(i,:))=1;
    end
    for i=1:size(MapIndex,1)
        M_ap(i,MapIndex(i,:))=1;
    end
    M_ap=char(M_ap);
    subs=[subs;M_ap];
end
%M_ap End
if (strcmp(BasisStr,'findf'))
    subs=[subs;find(f)];
end
subs=unique(subs,'rows');
end

function z=lambda_min(Q)
z=eig(Q);
z=min(z);
end