% Firstly, run the entry codes "Entry_for_PH_Due_to_LeftVentricularDiastolicDysfunction.m"
% Then, plot the curves of pressure over time

% Written by: Ziyin Dai, May 27, 2019, daiziyin@mail.dlut.edu.cn
% Corresponding author, Hong Tang, tanghong@dlut.edu.cn

Tall=0.7845; %The time of a cardiac cycle
step=0.0005; 
fs=1/step; 
samNum=ceil(Tall/step); 
t=zeros(5*samNum,1);
Plpap=zeros(5*samNum,1);
Plpad=zeros(5*samNum,1);
Plpv=zeros(5*samNum,1);
Prv=zeros(5*samNum,1);
Pla=zeros(5*samNum,1);

for num=1:samNum
    Plpap(num,:)=allP(10*samNum+num,20);
    Plpad(num,:)=allP(10*samNum+num,22);
    Plpv(num,:)=allP(10*samNum+num,24);
    Prv(num,:)=allP(10*samNum+num,18);
    Pla(num,:)=allP(10*samNum+num,25);
    num=num+1;
end

for num1=samNum+1:2*samNum
    Plpap(num1,:)=allP(200*samNum+num1,20);
    Plpad(num1,:)=allP(200*samNum+num1,22);
    Plpv(num1,:)=allP(200*samNum+num1,24);
    Prv(num1,:)=allP(200*samNum+num1,18);
     Pla(num1,:)=allP(200*samNum+num1,25);
    num1=num1+1;
end

for num2=2*samNum+1:3*samNum
    Plpap(num2,:)=allP(500*samNum+num2,20);
    Plpad(num2,:)=allP(500*samNum+num2,22);
    Plpv(num2,:)=allP(500*samNum+num2,24);
    Prv(num2,:)=allP(500*samNum+num2,18);
    Pla(num2,:)=allP(500*samNum+num2,25);
    num2=num2+1;
end


for num3=3*samNum+1:4*samNum
    Plpap(num3,:)=allP(700*samNum+num3,20);
    Plpad(num3,:)=allP(700*samNum+num3,22);
    Plpv(num3,:)=allP(700*samNum+num3,24);
    Prv(num3,:)=allP(700*samNum+num3,18);
    Pla(num3,:)=allP(700*samNum+num3,25);
    num3=num3+1;
end

for num4=4*samNum+1:5*samNum
    Plpap(num4,:)=allP(887*samNum+num4,20);
    Plpad(num4,:)=allP(887*samNum+num4,22);
    Plpv(num4,:)=allP(887*samNum+num4,24);
    Prv(num4,:)=allP(887*samNum+num4,18);
    Pla(num4,:)=allP(887*samNum+num4,25);
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
plot(t(v4)-t(v4(1)),Plpap(v4),'b-','linewidth',1.5) 
plot(t(v4)-t(v4(1)),Plpad(v4),'r-','linewidth',1.5)
plot(t(v4)-t(v4(1)),Plpv(v4),'b--','linewidth',1.5)
plot(t(v4)-t(v4(1)),Prv(v4),'k-','linewidth',1.5)
plot(t(v4)-t(v4(1)),Pla(v4),'r--','linewidth',1.5)

grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel(' Time in seconds','fontsize',10)
% legend('\itP\rm_l_p_a_p','\itP\rm_l_p_a_d','\itP\rm_l_p_v','\itP\rm_r_v','\itP\rm_l_a')
set(legend,'Orientation','horizontal');
axis([0  0.7845 0 80])
axis off
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)
