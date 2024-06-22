clc;clear;close all;

%% 使用说明
% 燕山大学YSU
% 机械原理B课程设计 
% 小型精刨机主切削运动机构的设计 机构2
% author:王湛然
% Email:wangzhanran@stumail.ysu.edu.cn

%% 基本数据条件
l1 = 92.4;
l6 = 350;
l3 = 757.4;
omegal1 = -2.4*pi;
H1 = 100;
G3 = 200;
G5 = 700;
J_s3 = 1.1;
l_s3 = l3/2;
g = 10;

phy = -(0:20:360)*pi/180;
si =[];
vi = [];
ai =[];
Mbi = [];
for i = 1:length(phy)
  [s,v,a,Mb] = sol_wzr(l1,l6,l3,omegal1,phy(i),G5,G3,g,l_s3,J_s3);
  si =[si;s];
  vi = [vi;v];
  ai =[ai;a];
  Mbi = [Mbi;Mb];
end
data = [si vi ai Mbi];
% 读取Excel文件中的数据
[~, ~, raw] = xlsread('all_data.xlsx');

% 指定填入数据的起始行、终止行、起始列、终止列
startRow = 2;
endRow = 19;
startCol = 7;
endCol = 10;

% 将数据填入指定的范围
for i = startRow:endRow
    for j = startCol:endCol
        raw{i, j} = data(i-startRow+1, j-startCol+1);
    end
end

% 将修改后的数据写回Excel文件
xlswrite('all_data.xlsx', raw);


function [s,v,a,Mb] = sol_wzr(l1,l6,l3,omegal1,phy1,G5,G3,g,l_s3,J_s3)
 %% 运动求解
    % b3移动量求解
    sb = sqrt((l1*cos(phy1))^2 + (l6+l1*sin(phy1))^2);
    % phy3转角求解
    phy3 = acos(l1*cos(phy1)/sb);
    % v23、omegal3求解
    sol1 = [cos(phy3) -sb*sin(phy3);...
            sin(phy3) sb*cos(phy3)]\(omegal1 * l1*[-sin(phy1);cos(phy1)]);
    omegal3 = sol1(2);
    v23 = sol1(1);
    % a23、alpha3求解
    sol2 = [cos(phy3) -sb*sin(phy3);...
            sin(phy3)  sb*cos(phy3)]\...
        (-[-omegal3*sin(phy3) -v23*sin(phy3)-sb*omegal3*cos(phy3);...
          omegal3*cos(phy3) v23*cos(phy3)-sb*omegal3*sin(phy3)]*sol1 ...
        - omegal1^2*l1*[cos(phy1);sin(phy1)]);
    alpha3 = sol2(2);
    % 刨刀 v a s求解
    vd = -omegal3*l3;
    v = vd*sin(phy3);
    a = -omegal3^2*l3*cos(phy3) - alpha3*l3*sin(phy3);
    s = l3*cos(phy3)+200;

%% 力分析
    % 刨刀架分析
    phy_cha = mod(phy1,-2*pi);
    if (phy_cha<-15.31*pi/180) && (phy_cha>-(180-15.31)*pi/180)
        F_r = 0;
    else 
        F_r = -4500;
    end
    R_34 = -F_r - G5*(-a)*(1e-3)/g;
    % 导杆机构力分析
    R_43 = - R_34;
    % 构件3基本参数
    a_s3x = -l_s3*(1e-3)*(omegal3^2*cos(phy3)+alpha3*sin(phy3));
    a_s3y = -l_s3*(1e-3)*(omegal3^2*sin(phy3)-alpha3*cos(phy3));
    P_I3x = -G3*a_s3x/g;
    P_I3y = -G3*a_s3y/g;
    P_I3 = sqrt(P_I3x^2+P_I3y^2);
    M_I3 = -J_s3*alpha3;
%     %构件3待求变量
%     syms R_23x R_23y R_63x R_63y 
%     %4个未知量，自然就是4个方程了啊
%     eq1 = P_I3x+R_23x+R_63x+R_43 == 0;
%     eq2 = P_I3y+R_23y+R_63y-G3 == 0;
%     eq3 = - R_43*l_s3*(1e-3)*sin(phy3) - R_23x*(sb-l_s3)*(1e-3)*sin(phy3) + R_63x*l_s3*(1e-3)*sin(phy3)...
%           + R_23y*(sb-l_s3)*(1e-3)*cos(phy3) - R_63y*l_s3*(1e-3)*cos(phy3) + M_I3== 0;% 对s3取矩,def逆时针+
%     eq4 = R_23x*cos(phy3)+R_23y*sin(phy3) == 0;
%     sol3 = solve(eq1,eq2,eq3,eq4);
%     R_23x = double(sol3.R_23x);
%     R_23y = double(sol3.R_23y);

    sol3 = [1 0 1 0;
            0 1 0 1;
            -(sb-l_s3)*(1e-3)*sin(phy3) (sb-l_s3)*(1e-3)*cos(phy3) l_s3*(1e-3)*sin(phy3) -l_s3*(1e-3)*cos(phy3);...
            cos(phy3) sin(phy3) 0 0]\[-P_I3x-R_43;
                                      -P_I3y+G3;
                                      R_43*l_s3*(1e-3)*sin(phy3)-M_I3;
                                      0];
    R_23x = sol3(1);
    R_23y = sol3(2);
    %构件1、2
    R_21x = -R_23x;
    R_21y = -R_23y;
    R_21 = sqrt(R_21x^2+R_21y^2);
    %对支架取矩,def逆时针+
    if phy_cha <=0 && phy_cha >=-pi/2
        Mb = -R_21x*l1*(1e-3)*abs(sin(phy1))-R_21y*l1*(1e-3)*abs(cos(phy1));
    end
    if phy_cha <-pi/2 && phy_cha >=-pi
        Mb = R_21x*l1*(1e-3)*abs(sin(phy1))-R_21y*l1*(1e-3)*abs(cos(phy1));
    end
    if phy_cha <-pi && phy_cha >=-3*pi/2
        Mb = -R_21x*l1*(1e-3)*abs(sin(phy1))-R_21y*l1*(1e-3)*abs(cos(phy1));
    end
    if phy_cha <-3*pi/2 && phy_cha >=-2*pi
        Mb = R_21x*l1*(1e-3)*abs(sin(phy1))-R_21y*l1*(1e-3)*abs(cos(phy1));
    end
end
