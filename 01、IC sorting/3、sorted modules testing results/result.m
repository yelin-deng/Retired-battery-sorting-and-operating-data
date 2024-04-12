%% 看看效果
% 同组电池的短期一致性（7天的电压和容量）和长期一致性（半年的电压和容量）；
% 与随机分配的对比。
clear;clc
path = 'C:\Users\11566\Desktop\毕业论文书写\第一章\materials\results and discussions';  % 读取数据
cd(path);    %把当前工作目录切换到指定文件夹
%%  计算分组的电压和容量的标准差
%加载数据
% 找最符合的，所有CAP和VOLT都引入。。
VOLT=load("VOLT.mat");
VOLT=VOLT.VOLT;
CAP=load("CAP.mat");
CAP=CAP.CAP;
%%
 c_k_group=load('c_k_group.mat');
c_k_group=c_k_group.c_k_group;
% c_k_group=c_k_group.group;

for i=1:length(c_k_group)
    onegroup=c_k_group{i};
    %电压
    vol=VOLT(onegroup,[2 end]);
    %容量
    cap=CAP(onegroup,[2 end]);
    ALL_CK=[vol cap];
%     ALL_CK=normalize(ALL_CK,'range');
    ALL_talk_CK(i,:)=std(ALL_CK);
    jiaquan(i)=length(onegroup);
end
   jiaquan1=jiaquan./sum(jiaquan);
   sum_in_CK=jiaquan1*ALL_talk_CK;
