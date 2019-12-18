% Firstly, run the entry codes "Entry_for_PH_Due_to_VentriculaSeptalDefect.m"
% Then, plot the curves of pressure-volume loops for the four heart chambers

% Written by: Ziyin Dai, May 27, 2019, daiziyin@mail.dlut.edu.cn
% Corresponding author, Hong Tang, tanghong@dlut.edu.cn

figure('unit','inch','position',[0.4 1 8 4])

%% P-V loop of left ventricle
v1=15/step+200:15/step+1800;  
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'r-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10)
axis([0 240 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=200/step-200:200/step+1800; 
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10)
axis([0 240 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=400/step-200:400/step+1800;  
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10)
axis([0 240 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=600/step-200:600/step+1800;  
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10)
axis([0 240 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=648/step-200:648/step+1800;  
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10) 
axis([0 240 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=688/step-200:688/step+1800;  
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10) 
axis([0 240 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=698/step-200:698/step+1800; 
subplot(2,2,1),hold on
plot(allVlv(v1),Plv(v1),'b-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨a£©Volume of LV (ml)','fontsize',10)
axis([0 240 0 150]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% P-V loop of left atrium
v1=15/step+200:15/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'r-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([20 200 0 60])  
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=200/step+200:200/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([20 200 0 60]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=300/step+200:300/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([20 200 0 60]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=400/step+200:400/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([20 200 0 60]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=600/step+200:600/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([20 200 0 60]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=680/step+200:680/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([20 200 0 60]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=698/step+200:698/step+1800;  
subplot(2,2,2),hold on
plot(allVla(v1),Pla(v1),'b-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨b£©Volume of LA (ml)','fontsize',10)
axis([20 200 0 60]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1) 

%% P-V loop of right ventricle
v1=15/step+200:15/step+1800;
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'r-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=150/step+200:150/step+1800; 
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=300/step+200:300/step+1800;  
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=400/step+200:400/step+1800;  
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=500/step+200:500/step+1800;  
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=600/step+200:600/step+1800; 
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=670/step+200:670/step+1800; 
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=680/step+200:680/step+1800; 
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=690/step+200:690/step+1800; 
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=698/step+200:698/step+1800;  
subplot(2,2,3),hold on
plot(allVrv(v1),Prv(v1),'b-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨c£©Volume of RV (ml)','fontsize',10)
axis([0 250 0 130]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

%% P-V loop of right atrium
v1=15/step+200:15/step+1800;  
subplot(2,2,4),hold on
plot(allVra(v1),Pra(v1),'r-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨d£©Volume of RA (ml)','fontsize',10)
axis([20 180 0 30]) 
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=200/step+200:200/step+1800; 
subplot(2,2,4),hold on
plot(allVra(v1),Pra(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨d£©Volume of RA (ml)','fontsize',10)
axis([20 180 0 30])  
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=400/step+200:400/step+1800;  
subplot(2,2,4),hold on
plot(allVra(v1),Pra(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨d£©Volume of RA (ml)','fontsize',10)
axis([20 180 0 30])
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=500/step+200:500/step+1800;  
subplot(2,2,4),hold on
plot(allVra(v1),Pra(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨d£©Volume of RA (ml)','fontsize',10)
axis([20 180 0 30])
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=600/step+200:600/step+1800; 
subplot(2,2,4),hold on
plot(allVra(v1),Pra(v1),'k-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨d£©Volume of RA (ml)','fontsize',10)
axis([20 180 0 30])
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)

v1=698/step+200:698/step+1800;  
subplot(2,2,4),hold on
plot(allVra(v1),Pra(v1),'b-','linewidth',1.5)
grid on
box on;
ylabel('BP in mmHg','fontsize',10)
xlabel('£¨d£©Volume of RA (ml)','fontsize',10)
axis([20 180 0 30])  
set(gca,'GridLineStyle','-.','GridColor','k', 'GridAlpha',1)
