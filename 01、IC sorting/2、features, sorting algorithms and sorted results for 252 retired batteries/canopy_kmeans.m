clear;clc
path = 'C:\Users\11566\Desktop\数据汇总\代码';  % 读取数据
cd(path);    %把当前工作目录切换到指定文件夹
%% 读取点

fea=load("4fea.mat");
fea=fea.IC;
fea_nored=normalize(fea,"range");
% IC_nored=IC_nored.*[0.1 0.1 0.5 0.1];
%% 读取分组数据
canopies_group=load('canopies_group.mat');
canopies_group=canopies_group.canopies_group;
%将canopies_group中的极小聚类删除
num_every_group=cellfun(@(x) length(x),canopies_group);
index_smallest=find(num_every_group<=5);
canopies_group_del=canopies_group(index_smallest);
canopies_group(index_smallest)=[];
%% 将canopies_group中的重复项从较少数量的簇中去除
full_group=cell(length(canopies_group),1);
full_group{1,1}=canopies_group{1,1};
for i=2:length(canopies_group)
 lastgroup=canopies_group(1:i-1);
 lastgroup= cellfun(@(x) x(:),lastgroup, 'UniformOutput', false);
 % 将剔除的聚类也放入其中，用于识别重复
 lastgroup_plus= cellfun(@(x) x(:),canopies_group_del, 'UniformOutput', false);
 lastgroup=cell2mat(lastgroup);
 lastgroup_plus=cell2mat(lastgroup_plus);
 lastgroup=[lastgroup;lastgroup_plus];
 nowgroup=canopies_group{i,1};
 %将last中有的从now中删除
 chongfu=intersect(lastgroup,nowgroup);
 full_group{i,1}=setdiff(nowgroup,chongfu);
end
% 可用的聚类
full_group;
% 剔除的聚类
canopies_group_del;
%% 开始进行kmeans
% 找出初始聚类中心
data0= cellfun(@(x) x(:),full_group, 'UniformOutput', false);
data=cell2mat(data0);
data=sort(data);
data=fea_nored(data,:);
k=length(full_group);
data_canopy=cellfun(@(x) fea_nored(x,:),full_group, 'UniformOutput', false); %获取点的数据
center=cellfun(@(x) mean(x),data_canopy, 'UniformOutput', false);      %计算质心
center=cell2mat(center);                                 %cell转化为矩阵
distance_center=pdist2(fea_nored,center);      % 计算每组点离质心的距离
[~,index_center]=min(distance_center);
center=fea_nored(index_center,:);
%%

idx=kmeans(data,5,"Start",center);
for id=1:5
 group{id,1}=find(idx==id);
end
c_value= cellfun(@(x) data(x,:), group , 'UniformOutput', false);
% end
% 计算group的数量
% jishu_group=cellfun(@(x) length(x),group, 'UniformOutput', false);
% jishu_group=cell2mat(jishu_group);
% sum(i,:)=sort(jishu_group);
% end
% z=mean(sum);
% sum2=sum-z;
% z2=std(sum2,0,2);
% [z4,z3]=min(z2);
%保存
c_k_group=group;
% 考虑剔除的电池，进行调整
% 得到剔除矩阵；矩阵中小于当前号码的数量；加上。
num_group=cellfun(@(x) length(x),group, 'UniformOutput', false);
num_group=cell2mat(num_group);
num_group=cumsum(num_group);% 用于后面的截取（索引）
num0= cellfun(@(x) x(:),group, 'UniformOutput', false);
num0=cell2mat(num0);% 用于后面的截取
% 剔除矩阵
num_del= cellfun(@(x) x(:),canopies_group_del, 'UniformOutput', false);
num_del=cell2mat(num_del);num_del=sort(num_del);
%加上一个矩阵
plus=zeros(length(fea_nored),1);
plus(num_del,:)=1;
plus=cumsum(plus);
plus(num_del,:)=[];
num1=num0;                         %num0为原始展开，num1为其增加数
for i_i=1:length(num1)
  work_num=num1(i_i);              % 找到应该增加的数
  num1(i_i)=work_num+plus(work_num);  % 执行增加
end
c_k_group{1}=num1(1:num_group(1));
for i=2:k
c_k_index=num_group(i-1)+1:num_group(i);
c_k_group{i}=num1(c_k_index);
end
%% 
save('c_k_group.mat','c_k_group');
save("canopies_group.mat",'full_group');
%% 将每个分组对应的IC-nored加上
c_k_group_value=c_k_group;     %canopy和kmeans联合方法
for i=1:length(c_k_group)

 work_group=c_k_group{i};
 
c_k_group_value{i}=fea_nored(work_group,:);

end
canopies_group_value=full_group;%canopy方法
for i=1:length(canopies_group)

 work_group=full_group{i};
 
canopies_group_value{i}=fea_nored(work_group,:);

end
kmean_group=load('kmean_group.mat',"group");
kmean_group=kmean_group.group;
 kmean_group_value=kmean_group;    %k means方法
for i=1:length(kmean_group)

 work_group=kmean_group{i};
 
kmean_group_value{i}=fea_nored(work_group,:);

end                                     
%% 计算SSE、MSE、最大距离
%SSE的计算：每个点与质心的距离的平方和
%MSE的计算：平均每个点与质心的距离
%最大距离的计算：计算每个点与质心的距离，取最大值
% work_cell=kmean_group_value;
% work_cell=canopies_group_value;
work_cell=c_k_group_value;
% work_cell=c_value;
for i=1:length(kmean_group_value)
    
    onegroup=work_cell{i};          %读取每一组
    meanpoint=mean(onegroup);       %计算质心
    distance_sum=pdist2(onegroup,meanpoint);     % 计算所有点到所有聚类中心的距离
    lg(i,1)=max(distance_sum);                   % 求最大得到半径
    sse0=distance_sum.^2;                        % 将距离取平方
    sse(i,1)=sum(sse0);                          % 求和得到sse
    mse(i,1)=sse(i,1)./length(onegroup);         % 求平均得到mse
    count(i,1)=length(onegroup);                 % 统计数量
end
count_weight=count./sum(count);
ave_lg=dot(lg,count_weight);
ave_sse=dot(sse,count_weight);
ave_mse=dot(mse,count_weight);









