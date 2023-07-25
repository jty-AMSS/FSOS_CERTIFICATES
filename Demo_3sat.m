File2opt=containers.Map('KeyType',  'char', 'ValueType', 'any');
File2opt('file_rwms_wcnf_L3_V70_C300_0.wcnf')=0
File2opt('file_rwms_wcnf_L3_V70_C300_1.wcnf')=1
File2opt('file_rwms_wcnf_L3_V70_C300_2.wcnf')=0
File2opt('file_rwms_wcnf_L3_V70_C300_3.wcnf')=0
File2opt('file_rwms_wcnf_L3_V70_C300_4.wcnf')=1
File2opt('file_rwms_wcnf_L3_V70_C300_5.wcnf')=0
File2opt('file_rwms_wcnf_L3_V70_C300_6.wcnf')=1
File2opt('file_rwms_wcnf_L3_V70_C300_7.wcnf')=0
File2opt('file_rwms_wcnf_L3_V70_C300_8.wcnf')=0
File2opt('file_rwms_wcnf_L3_V70_C300_9.wcnf')=0
Name=['file_rwms_wcnf_L3_V70_C300_0.wcnf'
'file_rwms_wcnf_L3_V70_C300_1.wcnf'
'file_rwms_wcnf_L3_V70_C300_2.wcnf'
'file_rwms_wcnf_L3_V70_C300_3.wcnf'
'file_rwms_wcnf_L3_V70_C300_4.wcnf'
'file_rwms_wcnf_L3_V70_C300_5.wcnf'
'file_rwms_wcnf_L3_V70_C300_6.wcnf'
'file_rwms_wcnf_L3_V70_C300_7.wcnf'
'file_rwms_wcnf_L3_V70_C300_8.wcnf'
'file_rwms_wcnf_L3_V70_C300_9.wcnf'
]
Table=[];
for i=1:10
    disp(i)
path=cd;
Path=[path, '\dir\' Name(i,:)];

tic;
[lb_p(i),~,~,~,~,Index]=HLM08(Path,inf,[],'p');
t_p=toc;
tic;


[f,w,~]=DICMS2function(Path);
opti=File2opt(Name(i,:));
d=2;
if opti~=0
f=1*f-opti;
end
f=f+0.5;
range=0:sum(w);
range=range-opti;
range=range(range>=0);
range=range(:);
range=range+0.5;
k=size(Index,1);
[Q,lb_my(i)]=FastLowerBound_Poly(f,d,range,k,length(f)*5);
t_lb=toc;
Table=[Table;i,t_p,lb_p(i),t_lb,lb_my(i)]
end