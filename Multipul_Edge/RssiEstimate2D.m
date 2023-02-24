%%%Rssi��λ�㷨 2D
clear;
Length=100;
Width=100;              %��ʼ������
Node_number=5;          %�۲�վ����������3��

% ���۲����ȷֲ���Ϊ�˺����ǽ���ƥ�䣬����ͬ�����ȷֲ�
% for i=1:Node_number
%     Node(i).x=Width*rand;
%     Node(i).y=Length*rand;%�۲�վλ�ó�ʼ��
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
Target.y=Length*rand;%Ŀ����ʵλ�ã����

%���۲�վ��Ŀ��̽��20�Σ�ȡƽ��ֵ��ΪRssiֵ
times_es=20;
% Z=[];%���۲�վ��Ŀ��̽��20��Rssiֵ
for i=1:Node_number
    for t=1:times_es
        [d]=Get_DIST(Node(i),Target);%�۲�վ��Ŀ�����ʵ����
        Rssi(i,t)=GetRssiValue(d,5);  %�õ�Rssi��ֵ
    end
end

ZZ=[];%����ʮ�ι۲��ƽ��ֵ
for i=1:Node_number
    ZZ(i)=sum(Rssi(i,:))/times_es;
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
    H=[H;2*(Node(i).x-Node(1).x),2*(Node(i).y-Node(1).y)];
    b=[b;Zd(1)^2-Zd(i)^2+Node(i).D-Node(1).D];
end
Estimate=((H'*H)\H')*b;%����Ŀ��λ��
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
line([Est_Target.x,Target.x],[Est_Target.y,Target.y],'Color','k');%ʵ��Ŀ�������Ŀ��֮�������
legend([h1,h2,h3],'�۲�վ','Ŀ��λ��','����λ��');
Error_Dist=Get_DIST(Est_Target,Target);
xlabel(['error=',num2str(Error_Dist),'m'])

%%%%%�Ӻ���
%������Ϊdʱ�����õõ�Rssi��ֵ
function value=GetRssiValue(d,Q)
    A=-42;n=2;%A,n�ڲ�ͬ��Ӳ��ϵͳȡֵ��һ��
    % Q=5;%�����������Rssi����ʱ�����ǳ���
    value=A-10*n*log10(d)+sqrt(Q)*randn;
end
%��Rssi��ֵ�������d
function d=GetDistByRssi(rssi)
    A=-42;n=2;%A,n�ڲ�ͬ��Ӳ��ϵͳȡֵ��һ��
    d=10^((A-rssi)/10/n);
end
function [dist]=Get_DIST(A,B)
    dist=sqrt((A.x-B.x)^2+(A.y-B.y)^2);
end