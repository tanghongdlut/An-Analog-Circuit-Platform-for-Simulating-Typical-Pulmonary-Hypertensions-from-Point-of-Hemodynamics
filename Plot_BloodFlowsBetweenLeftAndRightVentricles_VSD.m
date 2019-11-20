% Firstly, run the entry codes "Entry_for_PH_Due_to_VentriculaSeptalDefect.m"
% Then, plot the blood flow between ventricles

% Written by: Ziyin Dai, May 27, 2019, daiziyin@mail.dlut.edu.cn
% Corresponding author, Hong Tang, tanghong@dlut.edu.cn

Tall=0.7845; %The time of a cardiac cycle
step=0.0005; 
fs=1/step; 
samNum=ceil(Tall/step); 
t=zeros(5*samNum,1);
Qlv=zeros(5*samNum,1);
Qla=zeros(5*samNum,1); 
Qrv=zeros(5*samNum,1);
Qra=zeros(5*samNum,1); 
Qltor=zeros(5*samNum,1);

for num=1:samNum
    Qlv(num,:)=allQ(10*samNum+num,2);
    Qla(num,:)=allQ(10*samNum+num,1);
    Qrv(num,:)=allQ(10*samNum+num,26);
    Qra(num,:)=allQ(10*samNum+num,25);
    Qltor(num,:)=allQ(10*samNum+num,36);
    num=num+1;
end

for num1=samNum+1:2*samNum
    Qlv(num1,:)=allQ(300*samNum+num1,2);
    Qla(num1,:)=allQ(300*samNum+num1,1);
    Qrv(num1,:)=allQ(300*samNum+num1,26);
    Qra(num1,:)=allQ(300*samNum+num1,25);
    Qltor(num1,:)=allQ(300*samNum+num1,36);
    num1=num1+1;
end

for num2=2*samNum+1:3*samNum
    Qlv(num2,:)=allQ(500*samNum+num2,2);
    Qla(num2,:)=allQ(500*samNum+num2,1);
    Qrv(num2,:)=allQ(500*samNum+num2,26);
    Qra(num2,:)=allQ(500*samNum+num2,25);
    Qltor(num2,:)=allQ(500*samNum+num2,36);
    num2=num2+1;
end

for num3=3*samNum+1:4*samNum
    Qlv(num3,:)=allQ(700*samNum+num3,2);
    Qla(num3,:)=allQ(700*samNum+num3,1);
    Qrv(num3,:)=allQ(700*samNum+num3,26);
    Qra(num3,:)=allQ(700*samNum+num3,25);
    Qltor(num3,:)=allQ(700*samNum+num3,36);
    num3=num3+1;
end

for num4=4*samNum+1:5*samNum
    Qlv(num4,:)=allQ(887*samNum+num4,2);
    Qla(num4,:)=allQ(887*samNum+num4,1);
    Qrv(num4,:)=allQ(887*samNum+num4,26);
    Qra(num4,:)=allQ(887*samNum+num4,25);
    Qltor(num4,:)=allQ(887*samNum+num4,36);
    num4=num4+1;
end

figure('unit','inch','position',[0.4 1 1.6 2])
t=0:step:Tall*5;
 
%   v4=1:1569;
%   v4=1569+1:1569*2;
%   v4=2*1569+1:1569*3;
%   v4=3*1569+1:1569*4;
   v4=4*1569+1:1569*5;
hold on
plot(t(v4)-t(v4(1)),Qlv(v4),'b-','linewidth',1.5) 
plot(t(v4)-t(v4(1)),Qla(v4),'b--','linewidth',1.5)
plot(t(v4)-t(v4(1)),Qrv(v4),'k-','linewidth',1.5)
plot(t(v4)-t(v4(1)),Qra(v4),'k--','linewidth',1.5)
plot(t(v4)-t(v4(1)),Qltor(v4),'r-','linewidth',1.5)

grid on
box on;
ylabel('BF in ml/s','fontsize',10)
xlabel(' Time in seconds','fontsize',10)
% legend('\itQ\rm_l_v','\itQ\rm_l_a','\itQ\rm_r_v','\itQ\rm_r_a','\itQ\rm_l_t_o_r')
set(legend,'Orientation','horizontal');
axis([0  0.7845 -100 2000])
axis off
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)
