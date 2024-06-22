clc;clear;close all;

Fr = [];
for t = 0:0.01:30
    t = mod(t,5);
    theta = t*360/5;
    if theta > 15.31 && theta < (180-15.31)
        Fr= [Fr;0];
    else
        Fr= [Fr;-4500];
    end
end
Fr_t = [(0:0.01:30)',Fr];
% 指定CSV文件的文件名
filename = 'output.csv';
% 将矩阵数据写入CSV文件
writematrix(Fr_t, filename);
