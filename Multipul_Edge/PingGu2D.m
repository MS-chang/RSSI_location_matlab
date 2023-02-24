%%%Rssi定位算法 2D更改后进行绘图评估
clear;

% 与Q相关的性能分析
for i=1:100
    for j=1:10
        err(i,j)=Ping2D(i,20);
    end
    Err(i)=sum(err(i,:))/10;
end
figure(1);
plot(real(Err(4:100)));
title("Q");

% 与RSSI检测次数相关的性能分析
% for i=10:100
%     for j=1:10
%         err2(i,j)=Ping2D(4,i);
%     end
%     Err2(i)=sum(err2(i,:))/10;
% end
% figure(2);
% plot(real(Err2(10:100)));
% title("TimesRSSI");

function [Error_Dist] = Ping2D(Q,times_es)
    Length=100;
    Width=100;              %初始化场地
    Node_number=5;          %观测站个数，至少3个

    % for i=1:Node_number
    %     Node(i).x=Width*rand;
    %     Node(i).y=Length*rand;%观测站位置初始化
    %     Node(i).D=Node(i).x^2+Node(i).y^2;
    % end

    Node(1).x=0;
    Node(1).y=0;
    Node(1).D=Node(1).x^2+Node(1).y^2;
    Node(2).x=100;
    Node(2).y=0;
    Node(2).D=Node(2).x^2+Node(2).y^2;
    Node(3).x=50;
    Node(3).y=50;
    Node(3).D=Node(3).x^2+Node(3).y^2;
    Node(4).x=0;
    Node(4).y=100;
    Node(4).D=Node(4).x^2+Node(4).y^2;
    Node(5).x=100;
    Node(5).y=100;
    Node(5).D=Node(5).x^2+Node(5).y^2;

    Target.x=Width*rand;
    Target.y=Length*rand;%目标真实位置，随机

    %各观测站对目标探测20次，取平均值作为Rssi值
    % times_es=20;
    % Z=[];%各观测站对目标探测20次Rssi值
    for i=1:Node_number
        for t=1:times_es
            [d]=Get_DIST(Node(i),Target);%观测站与目标的真实距离
            Rssi(i,t)=GetRssiValue(d,Q);  %得到Rssi的值
        end
    end

    ZZ=[];%储存十次观测的平均值
    for i=1:Node_number
        ZZ(i)=sum(Rssi(i,:))/times_es;
    end

    %根据Rssi求观测距离
    Zd=[];%计算的距离
    for i=1:Node_number
        Zd(i)=GetDistByRssi(ZZ(i));
    end

    %根据观测距离用最小二乘法估计目标位置
    H=[];b=[];
    for i=2:Node_number
        %三角测边法公式
        H=[H;2*(Node(i).x-Node(1).x),2*(Node(i).y-Node(1).y)];  
        b=[b;Zd(1)^2-Zd(i)^2+Node(i).D-Node(1).D];
    end
    Estimate=((H'*H)\H')*b;%估计目标位置
    Est_Target.x=Estimate(1);Est_Target.y=Estimate(2);
    Error_Dist=Get_DIST(Est_Target,Target);
end



%%%%%子函数
%当距离为d时，采用得到Rssi的值
function value=GetRssiValue(d,Q)
    A=-42;n=2;%A,n在不同的硬件系统取值不一样
    %Q=5;%噪声方差，由于Rssi测量时噪声非常大
    value=A-10*n*log10(d)+sqrt(Q)*randn;
end
%由Rssi的值计算距离d
function d=GetDistByRssi(rssi)
    A=-42;n=2;%A,n在不同的硬件系统取值不一样
    d=10^((A-rssi)/10/n);
end
function [dist]=Get_DIST(A,B)
    dist=sqrt((A.x-B.x)^2+(A.y-B.y)^2);
end