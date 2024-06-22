
# 小型精刨机主切削运动机构的设计

## 1 声明

课程:YSU机械原理B课程设计(方案二)

author:王湛然(车卓22-1)

指导教师:潘登

date:2024/06/21

Email: wangzhanran@stumail.ysu.edu.cn

## 2 项目要求

本项目要求设计一个小型精刨机的主切削运动机构，要求主切削运动机构能使刨刀作往复直线运动。在工作行程时，刨刀速度要平稳，近可能接近等速，在空回行程时，刨刀快速退回，即要有急回作用。具有良好的传力性能，结构简单，动作轻便灵活。主要设计参数如下：

直线刨削行程:  H=400mm

工件刨削长度:   L=320MM

行程速比系数:   K=1.41

每分钟往复切削运动的次数 n=72次/min

切削阻力:  Fr=4500N

滑枕导路支承面离安装平面的距离需控制在1000mm之内。

## 3 机构图

机构图如下

<img src="fig\image.png" alt="alt text" width="50%" height="50%">

各项参数

<img src="fig\image-1.png" alt="alt text" width="20%" height="20%">

## 4 项目说明

1. 图解法分析机构运动学与动力学：具体文件在caxa文件夹下

2. 解析法分析机构运动学与动力学：分析原理在报告中
3. 解析法代码实现：项目源代码在Source Code下，通过了两种语言(MATLAB和python)的测试
4. py源码需要在Source Code\python中pip install -r requirements.txt
5. 代码最后编译的可执行文件在solve_by_wzr下，有MATLAB的同学可以运行solve_by_wzr\matlab目录下的JYBks_wzr_solve.p和JYBks_wzr_solve.exe，没有的可以运行solve_by_wzr\python\dist目录下的JYBks_wzr_solve.exe(由于github上传文件大小限制，放在了网盘里)
6. 百度网盘链接：<a href="https://pan.baidu.com/s/1Ic-XxBND0fuvD7sq03AMsA?pwd=mw60">https://pan.baidu.com/s/1Ic-XxBND0fuvD7sq03AMsA?pwd=mw60</a>
提取码：mw60
--来自百度网盘超级会员V5的分享
7. 简单的机构运动简图、装配体文件、不太成熟的motion分析在sw文件夹下
8. 整体设计报告为目录下的，仅供参考
9. 最后运动图像为

<img src="fig\image-2.png" alt="alt text" width="50%" height="50%">

10. 在此，感谢潘登老师的指导，帮助解决了设计代码时的疑惑，感谢组员李亮、秦子越等对代码设计时以身试数，帮助修改代码中的BUG，还有好兄弟车辆22-4的屌丝谢博鑫的帮助nuo
11. 感谢YSU提供的课设机会，可以对自己的理论知识进行检验和应用，锻炼了我这三脚猫的coding能力，也算之前成图的老底还能吃.
12. 对以后课设的同学的建议:多与老师沟通，不仅能收获机原的知识，还有人生指导nuo
