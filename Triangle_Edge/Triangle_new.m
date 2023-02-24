%%%%%%%%%%%%%%%%%该函数是改进后的三角质心定位
function [P] = TriangleNew(AA,BB,CC,dA,dB,dC)
    %A,B,C为三个选定的信标节点,节点坐标已知(为便于防真及验证,代码中采用的等边三角形)
    A = [AA.x,AA.y];
    B = [BB.x,BB.y];
    C = [CC.x,CC.y];
    
    % dA=;
    % dB=;%只需赋上相应的值即可(通过之前的那个程序将RSSI转化成距离即可)
    % dC=;
    
    %定义未知坐标x,y为符号变量
    syms x y;
    
    %距离方程,以信标节点为圆心,信标节点到未知节点的测量距离为半径作三个圆
    f1 = (A(1)-x)^2+(A(2)-y)^2-dA^2;
    f2 = (B(1)-x)^2+(B(2)-y)^2-dB^2;
    f3 = (C(1)-x)^2+(C(2)-y)^2-dC^2;
    
    %任两个方程联立,求任两圆交点
    
    s1 = solve(f1,f2); %求A,B两圆的交点
    s2 = solve(f2,f3); %求B,C两圆的交点
    s3 = solve(f1,f3); %求A,C两圆的交点
    
    %将结果(符号变量)转换为双精度数值
    
    x1 = double(abs(s1.x));
    y1 = double(abs(s1.y));
    x2 = double(abs(s2.x));
    y2 = double(abs(s2.y));
    x3 = double(abs(s3.x));
    y3 = double(abs(s3.y));
    %选择内侧的三个交点
    
    %两圆相交于两点,距第三个圆心近的为选定交点Pab,Pbc,Pac
    d1(1) = sqrt(abs((C(1)-x1(1))^2+(C(2)-y1(1))^2));
    d1(2) = sqrt(abs((C(1)-x1(2))^2+(C(2)-y1(2))^2));
    if d1(1) <= d1(2)
        Pab(1) = x1(1);
        Pab(2) = y1(1);
        dab=d1(1);
    else
        Pab(1) = x1(2);
        Pab(2) = y1(2);
        dab=d1(2);
    end
    
    d2(1) = sqrt(((A(1)-x2(1))^2+(A(2)-y2(1))^2));
    d2(2) = sqrt(((A(1)-x2(2))^2+(A(2)-y2(2))^2));
    
    if d2(1) <= d2(2)
        Pbc(1) = x2(1);
        Pbc(2) = y2(1);
        dbc=d2(1);
    else
        Pbc(1) = x2(2);
        Pbc(2) = y2(2);
        dbc=d2(2);
    end
    
    d3(1) = sqrt(((B(1)-x3(1))^2+(B(2)-y3(1))^2));
    d3(2) = sqrt(((B(1)-x3(2))^2+(B(2)-y3(2))^2));
    
    if d3(1) <= d3(2)
        Pac(1) = x3(1);
        Pac(2) = y3(1);
        dac=d3(1);
    else
        Pac(1) = x3(2);
        Pac(2) = y3(2);
        dac=d3(2);
    end

    %%%%进行系数的计算
    dA_m1=1/dab;%dA的倒数
    dB_AB_m1=1/dbc^(dab/dbc);
    dC_AC_m1=1/dac^(dab/dac);

    %求三个圆内侧三个交点Pab,Pbc,Pac的质心,即为未知节点P,完成定位
    P(1) = (Pab(1)*(dA_m1+dB_AB_m1)+Pbc(1)*(dB_AB_m1+dC_AC_m1)+Pac(1)*(dA_m1+dC_AC_m1))/(2*(dA_m1+dB_AB_m1+dC_AC_m1));
    P(2) = (Pab(2)*(dA_m1+dB_AB_m1)+Pbc(2)*(dB_AB_m1+dC_AC_m1)+Pac(2)*(dA_m1+dC_AC_m1))/(2*(dA_m1+dB_AB_m1+dC_AC_m1));;
end