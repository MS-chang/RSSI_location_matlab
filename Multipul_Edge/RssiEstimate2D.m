%%%Rssi定位算法 2D
clear;
Length=100;
Width=100;              %初始化场地
Node_number=5;          %观测站个数，至少3个

% 将观测点均匀分布，为了和三角进行匹配，采用同样均匀分布
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
times_es=20;
% Z=[];%各观测站对目标探测20次Rssi值
for i=1:Node_number
    for t=1:times_es
        [d]=Get_DIST(Node(i),Target);%观测站与目标的真实距离
        Rssi(i,t)=GetRssiValue(d,5);  %得到Rssi的值
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

%%%
figure;
hold on;grid on;box on;axis([0 Width 0 Length]);
for i=1:Node_number
    h1=plot(Node(i).x,Node(i).y,'ko','Markerface','g','MarkerSize',8);
    text(Node(i).x+2,Node(i).y,['Station',num2str(i)]);
end
h2=plot(Target.x,Target.y,'p','Markerface','r','MarkerSize',10);
h3=plot(Est_Target.x,Est_Target.y,'d','Markerface','r','MarkerSize',10);
line([Est_Target.x,Target.x],[Est_Target.y,Target.y],'Color','k');%实际目标与估计目标之间的连线
legend([h1,h2,h3],'观测站','目标位置','估计位置');
Error_Dist=Get_DIST(Est_Target,Target);
xlabel(['error=',num2str(Error_Dist),'m'])

%%%%%子函数
%当距离为d时，采用得到Rssi的值
function value=GetRssiValue(d,Q)
    A=-42;n=2;%A,n在不同的硬件系统取值不一样
    % Q=5;%噪声方差，由于Rssi测量时噪声非常大
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