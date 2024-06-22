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

disp('机械原理B课程设计之《小型精刨机主切削运动机构的设计》,机构2')
disp('author:王湛然')
disp('车卓22-1,第二组')
disp('指导教师:潘登')
disp("Email:wangzhanran@stumail.ysu.edu.cn")
disp("代码说明:x右为0度,涉及转动的,逆时针为正,笛卡尔坐标系")
disp('感谢组员李亮、秦子越、孔昱翔以身试数,帮助修改代码中的bug')
disp('推荐组长打开并行计算,嗯,已经打开了,嘿嘿嘿,之后会帮你关上的,保护好你的cpu诺')
disp('作为组长,总要为组员做点什么吧');
disp('Once you get the physics right, the rest is mathematics!')
%% 组长专属
ifzu = input('你是组长吗?,是为1,不是为0\n');
if ifzu == 1
ifzhouqi = input('周期绘图?,是为1,不是为0\n');

        if ifzhouqi == 1
            nzhouqi = input('几个周期?\n');
            v = [];
            a = [];
            s = [];
            Mb = [];
            for i = -nzhouqi*2:0.01:0
                phy1 = pi*i;
                [vi,ai,si,Mb_i] = sol_wzr(l1,l6,l3,omegal1,phy1,G5,G3,g,l_s3,J_s3);
                v = [v,vi];
                a = [a,ai];
                s = [s,si];
                Mb = [Mb,Mb_i];
            end  
            phy1_n = (-nzhouqi*2:0.01:0)*pi;
        end
        if ifzhouqi == 0
            v = [];
            a = [];
            s = [];
            Mb = [];
            for i = -2:0.001:0
                phy1 = pi*i;
                [vi,ai,si,Mb_i] = sol_wzr(l1,l6,l3,omegal1,phy1,G5,G3,g,l_s3,J_s3);
                v = [v,vi];
                a = [a,ai];
                s = [s,si];
                Mb = [Mb,Mb_i];
            end  
            phy1_n = (-2:0.001:0)*pi;
        end

%% 绘图
    % 设置背景颜色
    set(gcf, 'Color', [1 1 1]);  
    subplot(2,2,1)
    plot(phy1_n,v,'LineWidth',1.5,'Color','b')
    title("速度分析",'FontSize',14)
    xlabel('曲柄转角\theta_1/rad','FontSize',12)
    ylabel('刨刀速度v/(m/s)','FontSize',12)
    hold on;
    grid on;
    axis tight;
    set(gca, 'XDir', 'reverse');
    set(gca, 'FontSize', 12)
    legend('速度','Location','best')
    
    subplot(2,2,2)
    plot(phy1_n,a,'LineWidth',1.5,'Color','r')
    title("加速度分析",'FontSize',14)
    xlabel('曲柄转角\theta_1/rad','FontSize',12)
    ylabel('刨刀加速度a/(mm/s^2)','FontSize',12)
    hold on;
    grid on;
    axis tight;
    set(gca, 'XDir', 'reverse');
    set(gca,'FontSize',12)
    legend('加速度','Location','best')
    
    subplot(2,2,3)
    plot(phy1_n,s,'LineWidth',1.5,'Color','g')
    title("位移分析",'FontSize',14)
    xlabel('曲柄转角\theta_1/rad','FontSize',12)
    ylabel('刨刀位移s/mm','FontSize',12)
    hold on;
    grid on;
    axis tight;
    set(gca, 'XDir', 'reverse');
    set(gca,'FontSize',12)
    legend('位移','Location','best')
    
    subplot(2,2,4)
    plot(phy1_n,Mb,'LineWidth',1.5,'Color','k')
    title("平衡力矩分析",'FontSize',14)
    xlabel('曲柄转角\theta_1/rad','FontSize',12)
    ylabel('平衡力矩Mb/Nm','FontSize',12)
    hold on;
    grid on;
    axis tight;
    set(gca, 'XDir', 'reverse');
    set(gca,'FontSize',12)
    legend('平衡力矩','Location','best')
    disp('作为组长,总要为组员做点什么吧')
    disp('Once you get the physics right, the rest is mathematics!')
end
if ifzu == 0
    % 独属于你的phy
    phy_for_you = input('输入你的角度(角度制),逆时针为正\n');
    phy1_wzr = phy_for_you*pi/180;
    [v_my,a_my,s_my,Mb_my,omegal3_my,alpha3_my,R_34_my,P_I3_my,M_I3_my,R_21_my] = sol_wzr(l1,l6,l3,omegal1,phy1_wzr,G5,G3,g,l_s3,J_s3)
    disp('大部分的人的图解法与我的解析法应该都有误差,保持在2%以内即可!!!')
    disp('严格了一些,但是尽量满足吧')
    disp('Once you get the physics right, the rest is mathematics!')
    pause;
end

function [v,a,s,Mb,omegal3,alpha3,R_34,P_I3,M_I3,R_21] = sol_wzr(l1,l6,l3,omegal1,phy1,G5,G3,g,l_s3,J_s3)
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
