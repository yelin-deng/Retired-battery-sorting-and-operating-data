clear;clc
path = 'C:\Users\11566\Desktop\数据汇总\代码';  % 读取数据
cd(path);    %把当前工作目录切换到指定文件夹
%% 读取点

fea=load("4fea.mat");
fea=fea.IC;
fea_nored=normalize(fea,"range");

%% 计算所有点之间的距离
k2=1;
for i1=1:length(fea_nored)
    for j1=i1+1:length(fea_nored)
   distance_sum_points(k2,1)=norm(fea_nored(i1,:)-fea_nored(j1,:));
   k2=k2+1;
    end
end
sum_pingjun=mean(distance_sum_points);
sum_zhongwei=median(distance_sum_points);
sum_zuida=max(distance_sum_points);


T1=0.3;T2=0.28;
%% 聚类
% 初始化Canopies集合  
list=fea_nored;
list_index=1:length(list);
list_index=list_index';
canopies = [];  
canopies_group={};
% 随机选择一个点作为起始点  
start_point = randi(size(list,1));  
k1=1;
canopies(k1,:)  = list(start_point,:);  
canopies_group{k1} = start_point;  
list_index(start_point)=[];

%% 迭代处理数据点集合  
cishu=0;
while ~isempty(list_index)&&cishu<=10000%控制迭代次数
    cishu=cishu+1;
    % 随机从list中抽取一个
      i0=randi(size(list_index,1)); 
      i=list_index(i0);
      points=list(i,:);
    % 查找距离当前点最近的Canopy  
      [~,idx] = min(pdist2( points, canopies)); 
      distance=norm(points - canopies(idx,:) );
    % 判断
       if distance <= T2%归入该组，从list删除           
            canopies_group{idx,1} = union(canopies_group{idx}, i); 
            list_index=setdiff(list_index,i);
       elseif distance > T1  %新的canopy
             k1=k1+1;
            list_index=setdiff(list_index,i);
            canopies(k1,:) = points; %新坐标
            canopies_group{k1,1}=i;%新索引
       
       else
           canopies_group{idx,1} = union(canopies_group{idx}, i);  
       end
       
  


end  
%% 总结判断
   yifenzu= cellfun(@(x) x(:),canopies_group, 'UniformOutput', false);
   yifenzu=cell2mat(yifenzu);
   yifenzu=unique(yifenzu);
%% 生成聚类结果  
save('canopies_group.mat','canopies_group');

