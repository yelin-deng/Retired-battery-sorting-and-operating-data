%% 版本说明，与原版本相比，v的获取存在差距
% 原版本为对顺逆序电压分别计算ic，该版本在电压阶段对顺逆序耦合
% 并合理推测q的对应点
% 此外，由于后面的面积等特征所需的变量已被删除，所以将整块删除
%%
%选取第一次和第三次的特征
%IC峰和坐标，DV峰和坐标
%形成252行，4列的double数据
clear;clc
path = 'C:\Users\11566\Desktop\毕业论文书写\第一章\materials';  % 读取数据
cd(path);    %把当前工作目录切换到指定文件夹
%% 获取第一次循环的放电和第三次循环的数据
data=load("3cycle.mat");
data=data.data;
%% 开始循环
data0=data{1,1};%原始循环数据
for i=1:252
% 原始的时间、电流、电压
data1=data0{i};
t0=data1(1:3600,1);%s
I=data1(1:3600,2);%A
v0=data1(1:3600,3);%V
% 计算电池容量
t01=[0;t0];
t1=diff(t01);
q0=I.*t1./3600;%AH
q1=cumsum(q0);
%% 调节v0，使波动消除
%正序
v2_1=v0;
v0_a=unique(v2_1);
v0_s=mean(diff(v0_a));
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
% plot(v0,'y');hold on
% plot(v2_1-v0,'r'); hold on
% plot(v2_2-v0,'g'); hold on
% plot(v2-v0,'b');
%突然的想法
v3=(v2_2+v2_1)./2;% 顺逆序电压取平均
v3_a=unique(v3);% 找到不重复值，目的是避免ΔV=0
% 乱入计算平均步长
v3_s=mean(diff(v3_a));
% 乱入计算平均步长
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
q_a=q1(index_a_1);
ic3_aa=diff(q_a)./diff(v3_a);
v3_aa=v3_a(2:end);
%% 滤波——得先插值，再光滑smoothdata，
% 插值为与连续电压一一对应的，再光滑。
v4=min(v3_aa):0.001:max(v3_aa);
ic4 = interp1(v3_aa,ic3_aa,v4,'linear');%更优于三次样条插值
% plot(v3_aa,ic3_aa,'o',v4,ic4);
% 光滑
ic5=smoothdata(ic4,'lowess',6);
% plot(v4,ic4,':.',v4,ic5,'b');
% plot(v4,ic5);hold on
[~,ic_index]=max(ic5);
ic_y_high=ic5(ic_index);
ic_x_high=v4(ic_index);
% plot(ic_x_high,ic_y_high,'o');hold on
%% 根据观察，A峰在3.23-3.27V附近。代码的功能是寻找3.23V-3.27V的最大IC值和电压。
[~,work1]=min(abs(v4-3.23));
[~,work2]=min(abs(v4-3.27));
w=ic5(work1:work2);
[~,pIC] = max(w);
pIC=pIC+work1-1;
ic_y_low=ic5(pIC);%A峰IC值
ic_x_low=v4(pIC);%A峰电压
% plot(ic_x_low,ic_y_low,'*');hold on
mean_v=mean(v1(1:3200));
std_v=std(v1(1:3200));
%% 画图部分
% Q=qa;
% V(:,1)=v0;
%只能处理一下了，，
qz=0:0.5:q1(end);
qz=qz';
vz = interp1(q_a,v3_a,qz,'linear');%插值使数据长度一样
VZ(:,i)=vz';
plot(q_a,v3_a);hold on;
iz=interp1(q1,I,qz,'linear');

IC_skill(i,:)=[ic_x_low ic_y_low mean_v std_v];
end
% 保存IC特征集合
% save('IC_skill.mat','IC_skill');

