%%%%%%%%%%%%%%%%2022/4/6 Rssi定位算法
clc;
Length=100;
Width=100;
Hight=100;              %初始化场地
Node_number=5;          %观测站个数，至少3个
for i=1:Node_number
    Node(i).x=Width*rand;
    Node(i).y=Length*rand;
    Node(i).z=Hight*rand;
    %观测站位置初始化
    Node(i).D=Node(i).x^2+Node(i).y^2+Node(i).z^2;
end
Target.x=Width*rand;
Target.y=Length*rand;
Target.z=Hight*rand;   %目标真实位置，随机
%各观测站对目标探测10次，取平均值作为Rssi值
Z=[];%各观测站对目标探测10次Rssi值
for i=1:Node_number
    for t=1:20
        [d]=Get_DIST(Node(i),Target);%观测站与目标的真实距离
        Rssi(i,t)=GetRssiValue(d);  %得到Rssi的值
    end
end
ZZ=[];%储存十次观测的平均值
for i=1:Node_number
    ZZ(i)=sum(Rssi(i,:))/20;
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
    H=[H;2*(Node(i).x-Node(1).x),2*(Node(i).y-Node(1).y),2*(Node(i).z-Node(1).z)];  
    b=[b;Zd(1)^2-Zd(i)^2+Node(i).D-Node(1).D];
end
Estimate=inv(H'*H)*H'*b;%估计目标位置
Est_Target.x=Estimate(1);Est_Target.y=Estimate(2);Est_Target.z=Estimate(3);
%%%%%%%%%%%%
figure
for i=1:Node_number
    h1=plot3(Node(i).x,Node(i).y,Node(i).z,'ko','Markerface','g','MarkerSize',8);
    text(Node(i).x+2,Node(i).y,Node(i).z,['Station',num2str(i)]);
    hold on;
end
h2=plot3(Target.x,Target.y,Target.z,'p','Markerface','r','MarkerSize',10);
h3=plot3(Est_Target.x,Est_Target.y,Est_Target.z,'d','Markerface','r','MarkerSize',10);
line([Est_Target.x,Target.x],[Est_Target.y,Target.y],[Est_Target.z,Target.z],'Color','k');%实际目标与估计目标之间的连线
legend([h1,h2,h3],'观测站','目标位置','估计位置');
Error_Dist=Get_DIST(Est_Target,Target);
xlabel(['error=',num2str(Error_Dist),'m'])
box on;grid on;axis([0 Width 0 Length 0 Hight]);
%%%%%子函数
%当距离为d时，采用得到Rssi的值
function value=GetRssiValue(d)
    A=-42;n=2;%A,n在不同的硬件系统取值不一样
    Q=5;%噪声方差，由于Rssi测量时噪声非常大
    value=A-10*n*log10(d)+sqrt(Q)*randn;
end
%由Rssi的值计算距离d
function d=GetDistByRssi(rssi)
    A=-42;n=2;%A,n在不同的硬件系统取值不一样
    d=10^((A-rssi)/10/n);
end
function [dist]=Get_DIST(A,B)
    dist=sqrt((A.x-B.x)^2+(A.y-B.y)^2+(A.z-B.z)^2);    
end