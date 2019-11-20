% Firstly, run the entry codes "Entry_for_NormalHumanCirculationSystem.m"
% Then, plot the curves of pressure-volume loops for the four heart chambers 

% Written by: Ziyin Dai, May 27, 2019, daiziyin@mail.dlut.edu.cn
% Corresponding author, Hong Tang, tanghong@dlut.edu.cn

figure('unit','inch','position',[0.4 1 8 4])

%% The P-V loop of left ventricle
v1=698/step-200:698/step+1800; 
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10)
axis([20 160 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% The P-V loop of left atrium
v1=698/step+200:698/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([40 120 0 30]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% The P-V loop of right ventricle
v1=698/step+200:698/step+1800; 
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([20 160 0 30]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% The P-V loop of right atrium 
v1=698/step+200:698/step+1800; 
subplot(2,2,4),hold on
plot(allVra(v1),Pra(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨d£©Volume of RA (ml)','fontsize',10)
axis([40 120 0 30]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)
