import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 定义基本数据条件
l1 = 92.4
l6 = 350
l3 = 757.4
omegal1 = -2.4 * np.pi
H1 = 100
G3 = 200
G5 = 700
J_s3 = 1.1
l_s3 = l3 / 2
g = 10

def sol_wzr(l1, l6, l3, omegal1, phy1, G5, G3, g, l_s3, J_s3):
    phy1 = np.array(phy1)  # Ensure phy1 is numpy array for element-wise operations
    
    # 运动求解
    sb = np.sqrt((l1 * np.cos(phy1)) ** 2 + (l6 + l1 * np.sin(phy1)) ** 2)
    phy3 = np.arccos(l1 * np.cos(phy1) / sb)
    
    sol1 = np.linalg.solve(np.array([[np.cos(phy3), -sb * np.sin(phy3)],
                                     [np.sin(phy3), sb * np.cos(phy3)]]),
                           omegal1 * l1 * np.array([-np.sin(phy1), np.cos(phy1)]))
    omegal3 = sol1[1]
    v23 = sol1[0]
    
    sol2 = np.linalg.solve(np.array([[np.cos(phy3), -sb * np.sin(phy3)],
                                     [np.sin(phy3), sb * np.cos(phy3)]]),
                           -np.dot(np.array([[-omegal3 * np.sin(phy3), -v23 * np.sin(phy3) - sb * omegal3 * np.cos(phy3)],
                                            [omegal3 * np.cos(phy3), v23 * np.cos(phy3) - sb * omegal3 * np.sin(phy3)]]),
                                  sol1) - omegal1 ** 2 * l1 * np.array([np.cos(phy1), np.sin(phy1)]))
    alpha3 = sol2[1]
    
    vd = -omegal3 * l3
    v = vd * np.sin(phy3)
    a = -omegal3 ** 2 * l3 * np.cos(phy3) - alpha3 * l3 * np.sin(phy3)
    s = l3 * np.cos(phy3) + 200  
    
    # 力分析
    phy_cha = np.mod(phy1, -2 * np.pi)
    F_r = np.where((phy_cha < -15.31 * np.pi / 180) & (phy_cha > -(180 - 15.31) * np.pi / 180), 0, -4500)
    R_34 = -F_r - G5 * (-a) * (1e-3) / g
    
    R_43 = -R_34
    a_s3x = -l_s3 * (1e-3) * (omegal3 ** 2 * np.cos(phy3) + alpha3 * np.sin(phy3))
    a_s3y = -l_s3 * (1e-3) * (omegal3 ** 2 * np.sin(phy3) - alpha3 * np.cos(phy3))
    P_I3x = -G3 * a_s3x / g
    P_I3y = -G3 * a_s3y / g
    P_I3 = np.sqrt(P_I3x ** 2 + P_I3y ** 2)
    M_I3 = -J_s3 * alpha3
    
    sol3 = np.linalg.solve(np.array([[1, 0, 1, 0],
                                     [0, 1, 0, 1],
                                     [-(sb - l_s3) * (1e-3) * np.sin(phy3), (sb - l_s3) * (1e-3) * np.cos(phy3),
                                      l_s3 * (1e-3) * np.sin(phy3), -l_s3 * (1e-3) * np.cos(phy3)],
                                    [np.cos(phy3), np.sin(phy3), 0, 0]]),
                           np.array([-P_I3x - R_43,
                                     -P_I3y + G3,
                                     R_43 * l_s3 * (1e-3) * np.sin(phy3) - M_I3,
                                     0]))
    R_23x = sol3[0]
    R_23y = sol3[1]
    
    R_21x = -R_23x
    R_21y = -R_23y
    R_21 = np.sqrt(R_21x ** 2 + R_21y ** 2)
    
    if np.logical_and(phy_cha <= 0, phy_cha >= -np.pi / 2):
        Mb = -R_21x * l1 * (1e-3) * np.abs(np.sin(phy1)) - R_21y * l1 * (1e-3) * np.abs(np.cos(phy1))
    elif np.logical_and(phy_cha < -np.pi / 2, phy_cha >= -np.pi):
        Mb = R_21x * l1 * (1e-3) * np.abs(np.sin(phy1)) - R_21y * l1 * (1e-3) * np.abs(np.cos(phy1))
    elif np.logical_and(phy_cha < -np.pi, phy_cha >= -3 * np.pi / 2):
        Mb = -R_21x * l1 * (1e-3) * np.abs(np.sin(phy1)) - R_21y * l1 * (1e-3) * np.abs(np.cos(phy1))
    elif np.logical_and(phy_cha < -3 * np.pi / 2, phy_cha >= -2 * np.pi):
        Mb = R_21x * l1 * (1e-3) * np.abs(np.sin(phy1)) - R_21y * l1 * (1e-3) * np.abs(np.cos(phy1))
    
    return v, a, s, Mb, omegal3, alpha3, R_34, P_I3, M_I3, R_21

if __name__ == "__main__":
    print("version:1.1")
    print("date:2024.6.21")
    print("'指导教师:潘登'")    
    print("机械原理B课程设计之《小型精刨机主切削运动机构的设计》，机构2")
    print("author:王湛然")
    print("Email:wangzhanran@stumail.ysu.edu.cn")
    print("'指导教师:潘登'")
    print("代码说明:x右为0度,涉及转动的,逆时针为正,笛卡尔坐标系")
    print("感谢组员李亮、秦子越以身试数,帮助修改代码中的bug")
    print("作为组长,总要为组员做点什么吧")
    print("Once you get the physics right, the rest is mathematics!")

    ifzu = int(input('你是组长吗?,是为1,不是为0\n'))
    
    if ifzu == 1:
        ifzhouqi = int(input('周期绘图?,是为1,不是为0\n'))
        
        if ifzhouqi == 1:
            nzhouqi = int(input('几个周期?\n'))
            v = np.array([])
            a = np.array([])
            s = np.array([])
            Mb = np.array([])
            phy1_n = np.arange(-nzhouqi * 2 * np.pi, 0, 0.01 * np.pi)
            for i in phy1_n:
                [vi, ai, si, Mb_i, _, _, _, _, _, _] = sol_wzr(l1, l6, l3, omegal1, i, G5, G3, g, l_s3, J_s3)
                v = np.append(v, vi)
                a = np.append(a, ai)
                s = np.append(s, si)
                Mb = np.append(Mb, Mb_i)

        elif ifzhouqi == 0:
            v = np.array([])
            a = np.array([])
            s = np.array([])
            Mb = np.array([])
            phy1_n = np.arange(-2 * np.pi, 0, 0.001 * np.pi)
            for i in phy1_n:
                [vi, ai, si, Mb_i, _, _, _, _, _, _] = sol_wzr(l1, l6, l3, omegal1, i, G5, G3, g, l_s3, J_s3)
                v = np.append(v, vi)
                a = np.append(a, ai)
                s = np.append(s, si)
                Mb = np.append(Mb, Mb_i)

        # 加载中文字体
        font_path = 'C:\\Windows\\Fonts\\msyh.ttc'   # 微软雅黑字体路径
        my_font = fm.FontProperties(fname=font_path)

        # 绘图
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))
        axs[0, 0].plot(phy1_n, v, 'b', linewidth=1.5)
        axs[0, 0].set_title('速度分析', fontsize=14, fontproperties=my_font)
        axs[0, 0].set_xlabel('曲柄转角θ1/rad', fontsize=12, fontproperties=my_font)
        axs[0, 0].set_ylabel('刨刀速度v/(m/s)', fontsize=12, fontproperties=my_font)
        axs[0, 0].grid(True)
        axs[0, 0].invert_xaxis()

        axs[0, 1].plot(phy1_n, a, 'r', linewidth=1.5)
        axs[0, 1].set_title('加速度分析', fontsize=14, fontproperties=my_font)
        axs[0, 1].set_xlabel('曲柄转角θ1/rad', fontsize=12, fontproperties=my_font)
        axs[0, 1].set_ylabel('刨刀加速度a/(mm/s^2)', fontsize=12, fontproperties=my_font)
        axs[0, 1].grid(True)
        axs[0, 1].invert_xaxis()

        axs[1, 0].plot(phy1_n, s, 'g', linewidth=1.5)
        axs[1, 0].set_title('位移分析', fontsize=14, fontproperties=my_font)
        axs[1, 0].set_xlabel('曲柄转角θ1/rad', fontsize=12, fontproperties=my_font)
        axs[1, 0].set_ylabel('刨刀位移s/mm', fontsize=12, fontproperties=my_font)
        axs[1, 0].grid(True)
        axs[1, 0].invert_xaxis()

        axs[1, 1].plot(phy1_n, Mb, 'k', linewidth=1.5)
        axs[1, 1].set_title('平衡力矩分析', fontsize=14, fontproperties=my_font)
        axs[1, 1].set_xlabel('曲柄转角θ1/rad', fontsize=12, fontproperties=my_font)
        axs[1, 1].set_ylabel('平衡力矩Mb/Nm', fontsize=12, fontproperties=my_font)
        axs[1, 1].grid(True)
        axs[1, 1].invert_xaxis()

        plt.tight_layout()
        plt.show()
        a = input('按任意键退出..........')

    elif ifzu == 0:
        phy_for_you = float(input('输入你的角度(角度制),逆时针为正\n'))
        phy1_wzr = phy_for_you * np.pi / 180
        result = sol_wzr(l1, l6, l3, omegal1, phy1_wzr, G5, G3, g, l_s3, J_s3)
        v_my, a_my, s_my, Mb_my, omegal3_my, alpha3_my, R_34_my, P_I3_my, M_I3_my, R_21_my = result
        print(f'刨刀速度: {v_my}')
        print(f'刨刀加速度: {a_my}')
        print(f'刨刀位移: {s_my}')
        print(f'平衡力矩: {Mb_my}')
        print(f'角速度: {omegal3_my}')
        print(f'角加速度: {alpha3_my}')
        print(f'力R_34: {R_34_my}')
        print(f'力P_I3: {P_I3_my}')
        print(f'力矩M_I3: {M_I3_my}')
        print(f'力R_21: {R_21_my}')
        print('大部分的人的图解法与我的解析法应该都有误差, 保持在2%以内即可!!!')
        print('严格了一些, 但是尽量满足吧')
        print('Once you get the physics right, the rest is mathematics!')
        a = input('按任意键退出..........')
