%最近邻网络构建
clear; close all; clc  
cd("C:\Users\11566\Desktop\毕业论文书写\第三章\materials_new\cluster 算法\");
load("252_v_ic_dv.mat")
%% 距离计算,填入矩阵
DIST=zeros(252,252);
% hang=1:3019;
% data=V(hang,:);
data=DV;
% 欧式距离计算
for i=1:251
    vec1=data(:,i);
    for j=i+1:252
     vec2=data(:,j);
     DIST(i,j)=norm(vec1-vec2);
     DIST(j,i)=DIST(i,j);
    end
end
DIST([112 116],:)=[];
DIST(:,[112 116])=[];
%% 噪音去除，每个点到其他点的平均值
cell_sort_num=1:252;
for j=1:1
   mean_dist=mean(DIST);
   mean_dist=normalize(mean_dist,"range");
   std_dist=std(mean_dist);
   mean_mean_dist=mean(mean_dist);
   %三西格玛原则
   x=find(mean_dist>mean_mean_dist+3*std_dist);
   %剩余的电池序号和距离矩阵
   cell_sort_num=setdiff(cell_sort_num,x);
   DIST(x,:)=[];DIST(:,x)=[];
end
%% 每计算一次最近邻，输出最近邻网络包含的电池，剩余电池
ner_matrix=DIST;
for i=1:length(DIST)
  work_vec=DIST(i,:);
  [value,site]=sort(work_vec);
  ner_matrix(i,:)=site;
end
    work_ner_sum=ner_matrix(:,2:5);%9/10近邻之间___3近邻
    cell_sum=unique(work_ner_sum);
    num_in_net=length(cell_sum);
    cell_out_net=setdiff(1:length(DIST),cell_sum);
%% 得到每个网络
% 最近邻网络图
 work_ner_sum;
 ner_matrix=0;
   for i=1:length(DIST)
     [site_i_hang,site_i_lie]=find(work_ner_sum==i);
     ner_matrix(i,site_i_hang)=1;
     ner_matrix(site_i_hang,i)=1;
   end
% 寻找最近邻个数
cell_net=cell(length(DIST),1);
for i=1:length(DIST)
   NET1=[];k=1;
   diedaicishu=10;
      NET_0=[NET1;i];
      NET_1=[NET_0 unique(work_ner_sum(NET_0,:))]';
    while  length(NET_1)~=length(NET_0)
      NET_0=unique(NET_1);
      NET_1=unique(cat(1,NET_0,unique(work_ner_sum(NET_0,:))));
      k=k+1;
    end
 cell_net{i,1}=NET_1;
end
all_net=1:length(DIST);
del_net=[];j=1;
for i=1:length(DIST)
 net_i=cell_net{i,1};
 if isempty(intersect(net_i,del_net))
    remain_net=setdiff(all_net,net_i);
    all_net=remain_net;
    del_net=[del_net;net_i];
    net_record{j,1}=net_i;j=j+1;
 end
 
end
%% 选择初始聚类中心__没用上
% 思想，随机选择一个，然后从与该电池距离最远的20%中随机选，
% 再累加两个距离矩阵，从距离最远的20%中随机选。
c1=randi(252);
% c2
[~,dist_vect_y]=sort(DIST(c1,:));
c2_origin=dist_vect_y(round(0.8*length(dist_vect_y)):end);
c2=c2_origin(randi(length(c2_origin)));
c2_site=find(c2_origin==c2);
% c3
[~,dist_vect_y]=sort(DIST(c1,:)+DIST(c2,:));
c3_origin=dist_vect_y(round(0.8*length(dist_vect_y)):end);
c3=c3_origin(randi(length(c3_origin)));
c3_site=find(c3_origin==c3);
% c4
[~,dist_vect_y]=sort(DIST(c1,:)+DIST(c2,:)+DIST(c3,:));
c4_origin=dist_vect_y(round(0.8*length(dist_vect_y)):end);
c4=c4_origin(randi(length(c4_origin)));
c4_site=find(c4_origin==c4);
% c5
[~,dist_vect_y]=sort(DIST(c1,:)+DIST(c2,:)+DIST(c3,:)+DIST(c4,:));
c5_origin=dist_vect_y(round(0.8*length(dist_vect_y)):end);
c5=c5_origin(randi(length(c5_origin)));
c5_site=find(c5_origin==c5);
%% 将距离矩阵转化为坐标，方便后续使用matlab自带模糊C函数_输入所有点坐标和初始聚类中心
distance_matrix=DIST;
    [row, col] = size(DIST);  
    gram_matrix = zeros(row, col);  
    for i = 1:row  
        for j = 1:col  
            gram_matrix(i,j) = 0.5 * ((distance_matrix(1,j)^2) + (distance_matrix(i,1)^2) - (distance_matrix(i,j)^2));  
        end  
    end  
eps=1e-5;
% 计算特征值和特征向量  
[V, D] = eig(gram_matrix);  
% 初始化 select_feature_values 和 feature_vectors  
select_feature_values = diag(D); 
feature_vectors = V;  
% 删除特征值为0的特征向量  
while any(abs(select_feature_values) <eps)  
    [~, idx] = min(abs(select_feature_values));  
    select_feature_values(idx) = [];  
    feature_vectors(:, idx) = [];  
    D(:, idx)= [];  
    D(idx,: )= [];  
end  
% 计算最终坐标矩阵，通过将 feature_vectors 与 eye_matrix 的平方根相乘得到  
coordinates = feature_vectors * sqrt(D);
% 验证一下
coordinates=coordinates';
for i=1:row 
    vec1=coordinates(:,i);
    for j=1:col  
     vec2=coordinates(:,j);
     DIST_yz(i,j)=norm(vec1-vec2);
     DIST_yz(j,i)=DIST_yz(i,j);
    end
end
diff_matrix=DIST_yz-DIST;
error_corrd=max(diff_matrix(:));
coordinates=coordinates';
%% 已知坐标，开始聚类
% c1;c2;c3;c4;c5;
% Nc=coordinates([c1 c2 c3 c4 c5],:);
options = [1.1 NaN 1e-7 NaN];
[centers,U] = fcm(coordinates,5,options);
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);
index4 = find(U(4,:) == maxU);
index5 = find(U(5,:) == maxU);
%% 得到最终聚类序号
cell_sort_num;
group=cell(5,1);
group{1,1}=cell_sort_num(index1);
group{2,1}=cell_sort_num(index2);
group{3,1}=cell_sort_num(index3);
group{4,1}=cell_sort_num(index4);
group{5,1}=cell_sort_num(index5);
%% 保存
cell_del_num=setdiff(1:252,cell_sort_num);
 save('group_dv.mat','group','cell_del_num');