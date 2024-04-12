%获取252电池的v\IC\DV
clear
clc
close all;
cd("C:\Users\11566\Desktop\毕业论文书写\第三章\materials_new\base data\");
load('3cycle.mat','cycle_3');
data=cycle_3(1);
data_cell=data{1};
%% 计算出容量，以及每AH容量的位置
example=data_cell{1};
time_seq=example(:,1);
i_seq=example(:,2);
q=i_seq(1:end).*[1;diff(time_seq)]./3600;
Q=cumsum(q);
Q_max=abs(floor(Q(end)));
Q_site=zeros(length(Q_max),1);
for ii=1:Q_max
[~,Q_site(ii,1)]=min(abs(abs(Q)-ii));
end
%% 提取电压序列
for i=1:252
    work=data_cell{i};
    v0=work(:,3); 
%正序
  v2_1=v0;
  for i0=1:length(v0)
  v1=v0(1:i0);
  v2_1(i0,1)=max(v1);%处理后的电压序列
  end
%反序
  v2_2=v0;
  for i0=1:length(v0)
  v1=v0(end+1-i0:end);
  v2_2(length(v0)+1-i0,1)=min(v1);%处理后的电压序列
  end
  v2=(v2_1+v2_2)./2;
V(:,i)=v2;
end
%% 提取不重复电压序列的应有差值
for i=1:252
    v3=V(:,i);
    v3_a=unique(v3);% 找到不重复值，目的是避免ΔV=0
% 对不重复值插值
index_a=zeros(1,length(v3_a));
for i_a=1:length(v3_a)
    v3_aaa=v3_a(i_a);
    index_a(i_a)=find(v3==v3_aaa,1);
end% 首先找到顺序的
%之后找到平均位置，相当于顺逆序取平均
diff_index_a=diff(index_a);
index_a_1=index_a;
for ii_a=1:length(diff_index_a)-1
      index_a_1(ii_a)=floor((index_a(ii_a+1)+index_a(ii_a))./2);
end
%通过索引找到对应的Q
q_a=Q(index_a_1);
ic3_aa=diff(q_a)./diff(v3_a);
v3_aa=v3_a(2:end);
% 插值为与连续电压一一对应的，再光滑。
v4=min(v3_aa):0.001:max(v3_aa);
ic4 = interp1(v3_aa,ic3_aa,v4,'linear');
IC{i,1}=[v4;ic4];
%求DV
dv_site=Q_site;  % 每Ah容量的节点
dv_del=[];       % 删除  
for j=2:length(dv_site)
   qv1=v3(dv_site(j-1));
   qv2=v3(dv_site(j));
   if qv1==qv2
    dv_del=[dv_del j];
   end
end
dv_site(dv_del)=[];
dv0=diff([v3(1);v3(dv_site)])./diff([Q(1);Q(dv_site)]);
dvx=[Q(1);Q(dv_site(1:end-1))];
dv1 = interp1(dvx,dv0,1:Q_max-1,'linear');
dv1(isnan(dv1))=0;
DV(:,i)=dv1;
end
%% 保存电压、IC、DV
save('252_v_ic_dv.mat','V','IC','DV');









