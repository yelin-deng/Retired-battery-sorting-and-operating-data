%最近邻网络构建
clear; close all; clc  
cd("C:\Users\11566\Desktop\毕业论文书写\第三章\materials_new\cluster 算法\");
load("252_v_ic_dv.mat")
%% 距离计算,填入矩阵
DIST=zeros(252,252);
data=IC;
data=cellfun(@(x) x(2,:),data,'UniformOutput',false);
% 先归一化电压和IC
v_ic=cat(2,data{:});
v_ic_length=cellfun(@(x) length(x),data,'UniformOutput',false);
v_ic_norm=normalize(v_ic','scale');
%封装回cell
k=1;
for i=1:252
    length_num=v_ic_length{i};
    data{i}=v_ic_norm(k:k+length_num-1,:)';
    k=k+length_num;
end
% DTW距离计算
for i=1:251
    vec1=data{i};
    for j=i+1:252
     vec2=data{j};
     DIST(i,j)=dtw(vec1,vec2);
     DIST(j,i)=DIST(i,j);
    end
end
DIST([112 116],:)=[];
DIST(:,[112 116])=[];
%% 噪音去除，每个点到其他点的平均值
cell_sort_num=1:252;
for j=1:1
   mean_dist=mean(DIST);
   std_dist=std(mean_dist);
   mean_mean_dist=mean(mean_dist);
   %三西格玛原则，正好去掉了四个
   x=find(mean_dist>mean_mean_dist+3*std_dist);
   %剩余的电池序号和距离矩阵
   cell_sort_num=setdiff(cell_sort_num,x);
   DIST(x,:)=[];DIST(:,x)=[];
end
%% 思想，随机选择一个，然后从与该电池距离最远的5%中随机选，
% 再累加两个距离矩阵，从距离最远的5%中随机选。
for jj=224:224%1:length(DIST)

c1=jj;
% c2
[~,dist_vect_y]=sort(DIST(c1,:));
c2_origin=dist_vect_y(end:end);
c2=c2_origin(randi(length(c2_origin)));
c2_site=find(c2_origin==c2);
% c3
[~,dist_vect_y]=sort(DIST(c1,:)+DIST(c2,:));
c3_origin=dist_vect_y(end:end);
c3=c3_origin(randi(length(c3_origin)));
c3_site=find(c3_origin==c3);
% c4
[~,dist_vect_y]=sort(DIST(c1,:)+DIST(c2,:)+DIST(c3,:));
c4_origin=dist_vect_y(end:end);
c4=c4_origin(randi(length(c4_origin)));
c4_site=find(c4_origin==c4);
% c5
[~,dist_vect_y]=sort(DIST(c1,:)+DIST(c2,:)+DIST(c3,:)+DIST(c4,:));
c5_origin=dist_vect_y(end:end);
c5=c5_origin(randi(length(c5_origin)));
c5_site=find(c5_origin==c5);
%% 开始聚类
    % data: 数据矩阵，每行是一个数据点  
    % c: 聚类中心的数目  
    % m: 模糊参数，通常设置为2  
    % max_iter: 最大迭代次数  
    % tol: 收敛容差  
    % 初始化  
    n = length(DIST);  c=5;
    U = rand(c, n); 
    U = U ./ sum(U, 2); % 归一化
    m=1.1;
    centers_index=[c1 c2 c3 c4 c5];
    centers = data(centers_index, :);  
        % 计算隶属度矩阵  
        for j = 1:c  
            for i=1:n
            u_fenzi=DIST(centers_index(j),i);  
            u_fenmu1=DIST(centers_index(1),i);  
            u_fenmu2=DIST(centers_index(2),i);  
            u_fenmu3=DIST(centers_index(3),i);  
            u_fenmu4=DIST(centers_index(4),i);  
            u_fenmu5=DIST(centers_index(5),i); 
            mi=2/(m-1);
            U(j,i) = 1./((u_fenzi/u_fenmu1)^mi+(u_fenzi/u_fenmu2)^mi+(u_fenzi/u_fenmu3)^mi+(u_fenzi/u_fenmu4)^mi+(u_fenzi/u_fenmu5)^mi);  
            end  
        end
        %计算目标函数
        J=0;
        for j = 1:c  
            for i=1:n
                if isnan(U(j,i))
               continue
           
                end
        J=U(j,i).*DIST(centers_index(j),i)+J; 
            end  
        end
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);
index4 = find(U(4,:) == maxU);
index5 = find(U(5,:) == maxU);
matrix_reject(jj,2)=length(index5)+length(index4)+length(index3)+length(index2)+length(index1)-length(DIST);
matrix_reject(jj,1)=  J;
end
[z1,z2]=min(matrix_reject(:,1));
%% 得到最终聚类序号
cell_sort_num;
group=cell(5,1);
group{1,1}=cell_sort_num(index1);
group{2,1}=cell_sort_num(index2);
group{3,1}=cell_sort_num(index3);
group{4,1}=cell_sort_num(index5);
group{5,1}=cell_sort_num(index4);
% 
for i=1:4
    work=group{i};
    for j=i+1:5
       work2=group{j};
       jiaoji=intersect(work,work2);
        if ~isempty(jiaoji)
           for k=1:length(jiaoji)
       work(work==jiaoji(k))=[];
           end
       else
       end
    end
    group_new{i,1}=work;
end
 group_new{5,1}=group{5,1};
cell_del_num=setdiff(1:252,cell_sort_num);
group=group_new;
%% 
save('group_ic.mat','group','cell_del_num');