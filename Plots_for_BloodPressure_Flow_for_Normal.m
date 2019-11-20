% Firstly, run the entry codes "Entry_for_NormalHumanCirculationSystem.m"
% Then, plot the curves of blood pressure and flow with respect to time

% Written by: Ziyin Dai, May 27, 2019, daiziyin@mail.dlut.edu.cn
% Corresponding author, Hong Tang, tanghong@dlut.edu.cn

figure('unit','inch','position',[0.4 1 8 5])

%% Systemic blood pressures; 
v4=450/step-1069:450/step+2069;  
hold on
subplot(2,2,1),hold on
plot(t(v4)-t(v4(1)),Plv(v4),'k-','linewidth',1.5) 
plot(t(v4)-t(v4(1)),Paop(v4),'k-.','linewidth',1.5) 
plot(t(v4)-t(v4(1)),Pvc(v4),'k--','linewidth',1.5) 
plot(t(v4)-t(v4(1)),Pla(v4),'b-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel(' (a) Time in seconds','fontsize',10)
legend('\itP\rm_l_v','\itP\rm_a_o_p','\itP\rm_v_c','\itP\rm_l_a')
set(legend,'Orientation','horizontal');
axis([0 1.569 0 130])
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% Systemic blood flows 
v4=450/step-1069:450/step+2069;  
hold on
subplot(2,2,2),hold on
plot(t(v4)-t(v4(1)),allQ(v4,2),'k-','linewidth',1.5) 
plot(t(v4)-t(v4(1)),allQ(v4,13),'k-.','linewidth',1.5)  
plot(t(v4)-t(v4(1)),allQ(v4,24),'k--','linewidth',1.5) 
plot(t(v4)-t(v4(1)),allQ(v4,1),'b-','linewidth',1.5)
grid on
box on;
ylabel('Blood Flow in ml/s','fontsize',10)
xlabel(' (b) Time in seconds','fontsize',10)
legend('\itI\rm_l_v','\itI\rm_a_o_p','\itI\rm_v_c','\itI\rm_l_a')
set(legend,'Orientation','horizontal');
axis([0 1.569 0 900])
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% Pulmonary blood pressures
v4=450/step-1069:450/step+2069;  
hold on
subplot(2,2,3),hold on
plot(t(v4)-t(v4(1)),Prv(v4),'k-','linewidth',1.5)
plot(t(v4)-t(v4(1)),Pra(v4),'k-.','linewidth',1.5) 
plot(t(v4)-t(v4(1)),Plpap(v4),'k--','linewidth',1.5) 
plot(t(v4)-t(v4(1)),Plpv(v4),'b-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel(' (c)  Time in seconds','fontsize',10)
legend('\itP\rm_r_v','\itP\rm_r_a','\itP\rm_l_p_a_p','\itP\rm_l_p_v')
set(legend,'Orientation','horizontal');
axis([0 1.569 0 20])
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% pulmonary blood flows
v4=450/step-1069:450/step+2069;  
hold on
subplot(2,2,4),hold on
plot(t(v4)-t(v4(1)),allQ(v4,26),'k-','linewidth',1.5) 
plot(t(v4)-t(v4(1)),allQ(v4,25),'k-.','linewidth',1.5)  
plot(t(v4)-t(v4(1)),allQ(v4,30),'k--','linewidth',1.5) 
plot(t(v4)-t(v4(1)),allQ(v4,34),'b-','linewidth',1.5)
grid on
box on;
ylabel('Blood Flow in ml/s','fontsize',10)
xlabel(' (d) Time in seconds','fontsize',10)
legend('\itI\rm_r_v','\itI\rm_r_a','\itI\rm_l_p_a_p','\itI\rm_l_p_v')
set(legend,'Orientation','horizontal');
axis([0 1.569 0 1100])
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)
