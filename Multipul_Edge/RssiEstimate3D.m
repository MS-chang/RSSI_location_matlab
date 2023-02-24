%%%%%%%%%%%%%%%%2022/4/6 Rssi��λ�㷨
clc;
Length=100;
Width=100;
Hight=100;              %��ʼ������
Node_number=5;          %�۲�վ����������3��
for i=1:Node_number
    Node(i).x=Width*rand;
    Node(i).y=Length*rand;
    Node(i).z=Hight*rand;
    %�۲�վλ�ó�ʼ��
    Node(i).D=Node(i).x^2+Node(i).y^2+Node(i).z^2;
end
Target.x=Width*rand;
Target.y=Length*rand;
Target.z=Hight*rand;   %Ŀ����ʵλ�ã����
%���۲�վ��Ŀ��̽��10�Σ�ȡƽ��ֵ��ΪRssiֵ
Z=[];%���۲�վ��Ŀ��̽��10��Rssiֵ
for i=1:Node_number
    for t=1:20
        [d]=Get_DIST(Node(i),Target);%�۲�վ��Ŀ�����ʵ����
        Rssi(i,t)=GetRssiValue(d);  %�õ�Rssi��ֵ
    end
end
ZZ=[];%����ʮ�ι۲��ƽ��ֵ
for i=1:Node_number
    ZZ(i)=sum(Rssi(i,:))/20;
end
%����Rssi��۲����
Zd=[];%����ľ���
for i=1:Node_number
    Zd(i)=GetDistByRssi(ZZ(i));
end
%���ݹ۲��������С���˷�����Ŀ��λ��
H=[];b=[];
for i=2:Node_number
    %���ǲ�߷���ʽ
    H=[H;2*(Node(i).x-Node(1).x),2*(Node(i).y-Node(1).y),2*(Node(i).z-Node(1).z)];  
    b=[b;Zd(1)^2-Zd(i)^2+Node(i).D-Node(1).D];
end
Estimate=inv(H'*H)*H'*b;%����Ŀ��λ��
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
line([Est_Target.x,Target.x],[Est_Target.y,Target.y],[Est_Target.z,Target.z],'Color','k');%ʵ��Ŀ�������Ŀ��֮�������
legend([h1,h2,h3],'�۲�վ','Ŀ��λ��','����λ��');
Error_Dist=Get_DIST(Est_Target,Target);
xlabel(['error=',num2str(Error_Dist),'m'])
box on;grid on;axis([0 Width 0 Length 0 Hight]);
%%%%%�Ӻ���
%������Ϊdʱ�����õõ�Rssi��ֵ
function value=GetRssiValue(d)
    A=-42;n=2;%A,n�ڲ�ͬ��Ӳ��ϵͳȡֵ��һ��
    Q=5;%�����������Rssi����ʱ�����ǳ���
    value=A-10*n*log10(d)+sqrt(Q)*randn;
end
%��Rssi��ֵ�������d
function d=GetDistByRssi(rssi)
    A=-42;n=2;%A,n�ڲ�ͬ��Ӳ��ϵͳȡֵ��һ��
    d=10^((A-rssi)/10/n);
end
function [dist]=Get_DIST(A,B)
    dist=sqrt((A.x-B.x)^2+(A.y-B.y)^2+(A.z-B.z)^2);    
end