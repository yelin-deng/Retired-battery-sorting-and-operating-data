clear;clc
path = 'C:\Users\11566\Desktop\毕业论文书写\第一章\materials\聚类';  % 读取数据
cd(path);    %把当前工作目录切换到指定文件夹
%% 读取点
%% 读取点

IC_skill=load("IC_skill.mat");
IC_skill=IC_skill.IC_skill;
% 归一化
IC=IC_skill(:,[1 2 10 11]);
IC_nored=normalize(IC,"range");%使用归一化（0，1） 
%选择聚类中心
num_group=5;
sum_first=zeros(num_group,1);
sum_first(1,1)=randi([1 252]);
for  i=1:num_group
  nzero=find(sum_first~=0);
  num_nzero=length(nzero);
        if  num_nzero==num_group
           break
        end
  %% 计算距离聚类中心的距离
  k=1;plus=0;%弥补distance中的j1与真实j1的差距
  distance=zeros(252-num_nzero,1);
  for j=1:252
      if ismember(j,sum_first)
          plus=plus+1;
          continue
      end
  a=IC_nored(sum_first(1:num_nzero),:);
  %应对已经有两个或以上的聚类
  a=mean(a);
  b=IC_nored(j,:);
  distance(k,1)=norm(a-b);
  k=k+1;
  end
%% 转轮盘
  lunpan=sum(distance);
%   rng(232)
  changdu=rand()*lunpan;
  lunpan_new=0;%累加的初始长度
  j1=0;%初始值
  while lunpan_new<changdu
      j1=j1+1;
      lunpan_new=lunpan_new+distance(j1);
  end
  %选出下一个
sum_first(i+1,1)=j1+plus;
end
%% 将每个点分到最近的聚类中心
% load('compare.mat',"group",'center','sum_first');

distance_sum=9999*ones(252,num_group);
for i=1:252
%如果是聚类中心就撤出
    if ismember(i, sum_first)  
        continue
    end
    for j=1:num_group
distance_sum(i,j)=norm(IC_nored(i)-IC_nored(sum_first(j)));
    end
end
% 根据聚类矩阵分点
remain=setdiff(1:252,sum_first);
sum_first_cell=num2cell(sum_first);
 [~,min_all]=min(distance_sum,[],2);
to_remain=min_all(remain);
for i1=1:num_group
[to_group,~]=find(to_remain==i1);
sum_first_cell{i1}=[sum_first_cell{i1}; to_group];
end
%完成分组
%%  
% 重新计算每组的聚类中心（该组的平均值，大概率是一个 虚 点）
center=cell(num_group,1);
for j2=1:num_group
center_pre=sum_first_cell{j2};
center{j2,1}=mean(IC_nored(center_pre,:));
end
% 计算这些点到聚类中心的距离
for i2=1:252
a=IC_nored(i2,:);
  for j=1:num_group
        b=center{j};
      dist(i2,j)=norm(a-b); 
  end
end
% 将这些点重新分配
[~,min_all]=min(dist,[],2);
for i1=1:num_group
[to_group,~]=find(min_all==i1);
sum_first_cell{i1}=to_group;
end
% 再次计算聚类中心
center_new=cell(num_group,1);
for j2=1:num_group
center_pre=sum_first_cell{j2};
center_new{j2,1}=mean(IC_nored(center_pre,:));
end
% 重复，直到聚类中心不变
cishu=1;
while ~isequal(center_new,center)
center=center_new;
cishu=1+cishu;
for j2=1:num_group
center_pre=sum_first_cell{j2};
center{j2,1}=mean(IC_nored(center_pre,:));
end
% 计算这些点到聚类中心的距离
for i2=1:252
a=IC_nored(i2,:);
  for j=1:num_group
        b=center{j};
      dist(i2,j)=norm(a-b); 
  end
end
% 将这些点重新分配
[~,min_all]=min(dist,[],2);
for i1=1:num_group
[to_group,~]=find(min_all==i1);
sum_first_cell{i1}=to_group;
end
% 再次计算聚类中心
center_new=cell(num_group,1);
for j2=1:num_group
center_pre=sum_first_cell{j2};
center_new{j2,1}=mean(IC_nored(center_pre,:));
end

end
%% 
%加个大循环
%给出评价
distance_sum=cell(num_group,1);%距离汇总
for i=1:num_group
    center0=center{i};
    points_index=sum_first_cell{i};
    points=IC_nored(points_index,:);
    pre_distance=(center0-points)';
    distance_sum{i,1}=vecnorm(pre_distance);
    distance_ave0=vecnorm(pre_distance);
    distance_ave(i,1)=mean(distance_ave0);%簇内平均距离
    distance_de(i,1)=max(distance_ave0);%簇内半径
    k2=1;
for i1=1:length(points)
    for j1=i1+1:length(points)
   distance_sum_points(k2,1)=norm(points(i1,:)-points(j1,:));
   k2=k2+1;
    end
end
   distance_lgd(i,1)=max(distance_sum_points);
   points_num(i,1)=length(points);
    
end
%% 输出分类
group=sum_first_cell;
save('kmean_group.mat',"group");
%%

basegroup=group;








      