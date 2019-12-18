% This script is a simulation program for the normal human circulation system 
% based on an analog circuit platform.
% The circuit platform is shown in a document  "Figure_analog_circuit_platform.doc"

% The simulation time is 700s. A cardiac cycle is 0.7845s (heart rate is about 76.5 beat per minute). 
% The time step size in numerical solution is 0.0005s. 
% The total blood volume in the circulation system is 4711 ml.
% The sympathetic frequencies and vagal frequency is set as 0.5. 
% The initial blood volume of each capacitor, current of each inductor, initial values of capacitances, 
% inductances and resistances in the model given in Appendix A. 
% The authors assume that the time-varying parameters is constant within a cardiac cycle 
% and has an increment or reduction between adjacent cycles.

% Written by: Ziyin Dai, May 27, 2019, daiziyin@mail.dlut.edu.cn
% Corresponding author, Hong Tang, tanghong@dlut.edu.cn

clc; close all; clear all
tic    % to start a timer 
 
%%  Load the initial values of resistances
%   Rm     Ra   Rhaa  Rlna  Rlca  Raop  Rrula  Rrica  Rlica  Rlula  Rsap   Rrsv
R=[0.02   0.02   8     12    12   1.2    0.5    0.5    0.5    0.5    0.5   0.27 ...  
    0.25  0.25   0.25   0.2   0.04   0.03  0.01   0.05    0.05   0.06   0.06   0.07   0.07];
%  Rrijv  Rlijv  Rlsv   Rsv    Rvc    Rt    Rp    Rrpap   Rlpap   Rrpad  Rlpad  Rrpv   Rlpv   

%% Load the initial values of compliances    
%  Chaa   Clna   Clca   Caop  Crula  Crica  Clica   Clula  Csap   Crsv    
C=[0.7    0.7    0.7     0.8    3      2      3      2      5       9 ...
    9      9     9    10    15     1.5     1.5     9      9      15     15 ];
% Crijv  Clijv  Clsv  Csv   Cvc   Crpap   Clpap  Crpad  Clpad   Crpv    Clpv

%% Load the initial values of viscoelastic resistances
%  Rchaa  Rclna  Rclca  Rcaop   Rcrpap  Rclpap  Rcrpad   Rclpad   Rrpv    Rlpv
Rc=[0.01  0.01   0.01   0.01    0.005   0.005   0.005    0.005    0.005  0.005 ];
 
%% Load the initial values of inductances 
%  Laop   Lrpap  Llpap   
L=[0.001  0.001  0.001 ];

%---yinit:  1     2      3      4      5     6       7     8        9     10    
%----------Vlv   Vhaa   Vlna   Vlca  Vaop  Vrula   Vrica  Vlica   Vlula   Vsap   
%----------11     12     13    14     15     16      17    18      19      20
%----------Vrsv  Vrijv  Vlijv  Vlsv   Vsv    Vvc    Vra     Vrv    Vrpa    Vlpa  
%----------21     22    23
%----------Vrpv  Vlpv   Vla   

%% Load initial values at t=0s
% Include the initial conditions of the blood volumes in four chambers, vessels and the blood flows in systemic and pulmonary aorta flows.
load yinit_normal.mat       
yinit=yinit_normal; 

%% Adjustable part             
Tall=700;    % Simulation  duration, time in seconds
step=0.0005; % Simulation step size, time in seconds
fs=1/step;   % Frequency
samNum=ceil(Tall/step);  % Number of points required for simulation
Hr=76*ones(samNum+2,1);  % Heart rate

%% Normalized  Vagal and sympathetic  frequencies
FhrV=0.5*ones(samNum+2,1); % Fhrv is normalized vagal frequency
FhrS=0.5*ones(samNum+2,1); % Fhrs is normalized sympathetic frequency 
Fcon=0.5*ones(samNum+2,1); % Fcon is sympathetic efferent discharge frequency 
Fvaso=0.5*ones(samNum+2,1);% Fvaso is normalized sympathetic efferent frequency

%-------------%
beatNum=1;   % Number of heartbeats
HrT=zeros(1,2000);  % Time of each cardiac cycle
num=1;  % Record the number of simulation steps

allP=zeros(samNum+1,25);  % All blood pressure values at each moment
allV=zeros(samNum+1,25);  % All volume values at each moment
allD=zeros(samNum+1,17);  % Valve status at each moment (1: open; 0: close)
allQ=zeros(samNum+1,35);  % All blood flow values at each moment

LVSV=zeros(893,1); % Left ventricular stroke volume
RVSV=zeros(893,1); % Right ventricular stroke volume

Rtime=0;

for t=0:step:Tall
    
    %% Obtain the solution of the previous cycle equation  
    Vlv=yinit(1);     Vhaa=yinit(2);     Vlna=yinit(3);     Vlca=yinit(4); 
    Vaop=yinit(5);    Vrula=yinit(6);    Vrica=yinit(7);    Vlica=yinit(8);
    Vlula=yinit(9);   Q10=yinit(10);     Vsap=yinit(11);    Vrsv=yinit(12);   
    Vrijv=yinit(13);  Vlijv=yinit(14);   Vlsv=yinit(15);    Vsv=yinit(16);    
    Vvc=yinit(17);    Vra=yinit(18);     Vrv=yinit(19);     Vrpap=yinit(20);   
    Vlpap=yinit(21);  Q24=yinit(22);     Q25=yinit(23);     Vrpad=yinit(24);   
    Vlpad=yinit(25);  Vrpv=yinit(26);    Vlpv=yinit(27);    Vla=yinit(28); 
      
    allVlv(num)=Vlv; allVla(num)=Vla;
    allVrv(num)=Vrv; allVra(num)=Vra;
   
    %% The P-V relationship of four heart chambers
    ttemp=t-sum(HrT(1:beatNum)); % At the current moment of a new cardiac cycle
    Plv=mycallP_normal(Vlv,ttemp,Fcon(num),1); % Left ventricle
    Pla=mycallP_normal(Vla,ttemp,Fcon(num),3); % Left atrium
    Prv=mycallP_normal(Vrv,ttemp,Fcon(num),2); % Right ventricle
    Pra=mycallP_normal(Vra,ttemp,Fcon(num),4); % Right atrium 
    
    %% P-V relationships of linear vessels
    Phaa=Vhaa/C(1);     Plna=Vlna/C(2);      Plca=Vlca/C(3);
    Paop=Vaop/C(4);     Prula=Vrula/C(5);    Prica=Vrica/C(6);
    Plica=Vlica/C(7);   Plula=Vlula/C(8);    Prsv=Vrsv/C(10);
    Prijv=Vrijv/C(11);  Plijv=Vlijv/C(12);   Plsv=Vlsv/C(13);   
    Prpap=Vrpap/C(16);  Plpap=Vlpap/C(17);   Prpad=Vrpad/C(18);  
    Plpad=Vlpad/C(19);  Prpv=Vrpv/C(20);     Plpv=Vlpv/C(21);  
  
   %% Nonlinear P-V relationships for specified vessels
   Kc=1000; Vsap_min=210; N0=50;
   Kp1=0.03; Kp2=0.2; Taop=0.1;
   Psap_a=Kc*log10((Vsap-Vsap_min)/N0+1);  
   Psap_p=Kp1*exp(Taop*(Vsap-Vsap_min))+Kp2*(Vsap-Vsap_min).^2; 
   Psap=Fvaso(num)*Psap_a+(1-Fvaso(num))*Psap_p;  % Proximal arteries  
   
   Kv=40; Vsv_max=3500;
   Psv=-Kv*log10((Vsv_max/Vsv) -0.99); % Systemic veins

   N1=0; N2=-5; K1=0.15; K2=0.4; Vvc_min=50; Vvc_0=130;  
    if  Vvc>=130   % Vena cava
        Pvc=N1+K1*(Vvc-Vvc_0);
    else
        Pvc=N2+K2*exp(Vvc/Vvc_min);
    end
       
    allP(num,:)=[Plv Phaa Plna Plca Paop Prula Prica Plica Plula Psap Prsv Prijv Plijv Plsv Psv Pvc Pra Prv Prpap Plpap Prpad Plpad Prpv Plpv Pla];
    allV(num,:)=[Vlv Vhaa Vlna Vlca Vaop Vrula Vrica Vlica Vlula Vsap Vrsv Vrijv Vlijv Vlsv Vsv Vvc Vra Vrv Vrpap Vlpap Vrpad Vlpad Vrpv Vlpv Vla];
      
    Kr=0.04; Vsap_max=250;
    R(11)=Kr*exp(4*Fvaso(num))+Kr*(Vsap_max/Vsap).^2;  % Rsap
    
    KR=0.001; Vvc_max=350; R0=0.025;
    R(17)=KR*(Vvc_max/Vvc).^2+R0;  % Rvc
  
    %% The values of all blood flows
    Q1=(Pla-Plv)/R(1);    
    Q3=(Plv-Phaa)/R(2);        Q4=(Plv-Plna)/R(2);
    Q5=(Plv-Plca)/R(2);        Q6=(Plv-Paop)/R(2);
    Q71=(Phaa-Prula)/R(3);     Q72=(Phaa-Prica)/R(3);
    Q8=(Plna-Plica)/R(4);      Q9=(Plca-Plula)/R(5);    
    Q11=(Prula-Prsv)/R(7);     Q12=(Prica-Prijv)/R(8);
    Q13=(Plica-Plijv)/R(9);    Q14=(Plula-Plsv)/R(10);  
    Q15=(Psap-Psv)/R(11);
    Q16=(Prsv-Pvc)/R(12);      Q17=(Prijv-Pvc)/R(13);
    Q18=(Plijv-Pvc)/R(14);     Q19=(Plsv-Pvc)/R(15);
    Q20=(Psv-Pvc)/R(16);  
    Q21=(Pvc-Pra)/R(17);       Q22=(Pra-Prv)/R(18);
    Q240=(Prv-Prpap)/R(19);    Q250=(Prv-Plpap)/R(19);
    Q24=(Prpap-Prpad)/R(20);   Q25=(Plpap-Plpad)/R(21);
    Q26=(Prpad-Prpv)/R(22);    Q27=(Plpad-Plpv)/R(23);
    Q28=(Prpv-Pla)/R(24);      Q29=(Plpv-Pla)/R(25);
       
    %% Valves  model  (1: open; 0: close)
    Dm=Pla>Plv;  
    D1=Plv>Phaa;  D2=Plv>Plna;
    D3=Plv>Plca;  D4=Plv>Paop;  
    Da=D1|D2|D3|D4;
    D51=Phaa>Prula; D52=Phaa>Prica; 
    D53=Prsv>Pvc; D54=Prijv>Pvc; 
    
    D6=Plijv>Pvc;  D7=Plsv>Pvc;   D8=Psv>Pvc;  
    Dt=Pra>Prv; 
    D9=Prv>Prpap;   D10=Prv>Plpap; 
    Dp=D9|D10;
   
    Q2=D1*Q3+D2*Q4+D3*Q5+D4*Q6; 
    Q7=D51*Q71+D52*Q72;  
    Q77=D53*Q16+D54*Q17;
    Q23=D9*Q240+D10*Q250;
    Q30=Q28+Q29;
    
   %% Solution of differential equations
   % Convert differential equations into difference equations,
   % V(t+step)-V(t)=step*Q(t), Q(t+step)-Q(t)=step*P(t)/L,
    
    %% Recursive calculation 
    Eright(1)=Dm*Q1-Da*Q2;
    Eright(2)=D1*Q3-D51*Q71-D52*Q72;
    Eright(3)=D2*Q4-Q8;
    Eright(4)=D3*Q5-Q9;
    Eright(5)=D4*Q6-Q10;
    Eright(6)=D51*Q71-Q11;
    Eright(7)=D52*Q72-Q12;
    Eright(8)=Q8-Q13;
    Eright(9)=Q9-Q14;
    Eright(10)=(Paop-Q10*R(6)-Psap)/L(1); 
    Eright(11)=Q10-Q15;
    Eright(12)=Q11-D53*Q16;
    Eright(13)=Q12-D54*Q17;
    Eright(14)=Q13-D6*Q18;
    Eright(15)=Q14-D7*Q19;
    Eright(16)=Q15-D8*Q20;
    Eright(17)=D53*Q16+D54*Q17+D6*Q18+D7*Q19+D8*Q20-Q21;
    Eright(18)=Q21-Dt*Q22;
    Eright(19)=Dt*Q22-Dp*Q23;
    Eright(20)=D9*Q240-Q24;
    Eright(21)=D10*Q250-Q25;
    
    Eright(22)=(Prpap-Q24*R(20)-Prpad)/L(2); 
    Eright(23)=(Plpap-Q25*R(21)-Plpad)/L(3); 
    Eright(24)=Q24-Q26;
    Eright(25)=Q25-Q27;
    Eright(26)=Q26-Q28;
    Eright(27)=Q27-Q29;
    Eright(28)=Q28+Q29-Dm*Q1; 
      
    allQ1(num)=Q2*Da;
    allQ2(num)=Q1*Dm;
    
    Q1=Q1*Dm;
    Q2=Q2*Da;
    Q3=Q3*D1;
    Q4=Q4*D2;
    Q5=Q5*D3;
    Q6=Q6*D4;
    Q71=Q71*D51;
    Q72=Q72*D52;
    Q16=Q16*D53;
    Q17=Q17*D54;
    Q18=Q18*D6;
    Q19=Q19*D7;
    Q20=Q20*D8;
    Q22=Q22*Dt;
    Q23=Q23*Dp;
    Q240=Q240*D9;
    Q250=Q250*D10;

    allD(num,:)=[Da Dm Dp Dt D1 D2 D3 D4 D51 D52 D53 D54 D6 D7 D8 D9 D10];
    allQ(num,:)=[Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q71 Q72 Q77 Q8 Q9 Q10 Q11 Q12 Q13 Q14 Q15 Q16 Q17 Q18 Q19 Q20 Q21 Q22 Q23 Q240 Q250 Q24 Q25 Q26 Q27 Q28 Q29 Q30];

    yinit=yinit+step*Eright;
    num=num+1;
     
    %% Heart rate control
    if ttemp>60/Hr(num-1)      
        Rtime=[Rtime sum(HrT)];
        beatNum=beatNum+1;
        HrT(beatNum-1)=ttemp;

       %% Left ventricular stroke volume
        LVSV(beatNum-1,1)=max(allV(num-floor(HrT(beatNum-1)/step):num-1,1))-min(allV(num-floor(HrT(beatNum-1)/step):num-1,1));
       %% Right ventricular stroke volume
        RVSV(beatNum-1,1)=max(allV(num-floor(HrT(beatNum-1)/step):num-1,18))-min(allV(num-floor(HrT(beatNum-1)/step):num-1,18));
        
       %% Calculate the mean proximal right pulmonary artery pressure (mPAP=1/3sPAP+2/3dPAP)
        sPAP(beatNum-1,1)=max(allP(num-floor(HrT(beatNum-1)/step):num-1,19)); % Systolic blood pressure of the proximal right pulmonary artery 
        dPAP(beatNum-1,1)=min(allP(num-floor(HrT(beatNum-1)/step):num-1,19)); % Diastolic blood pressure of the proximal right pulmonary artery
        mPAP(beatNum-1,1)=(1/3)*sPAP(beatNum-1,1)+(2/3)*dPAP(beatNum-1,1); % Mean proximal right pulmonary artery pressure
  
    end  
    
    %% Heart rate control
    Hr(num)=35+FhrS(num)*140-40*FhrS(num).^2-32*FhrV(num)+10*FhrV(num).^2-20*FhrV(num).*FhrS(num);  
 
end

toc 

HrT=HrT(find(HrT~=0));
HrT=[HrT ttemp];

Plv=allP(:,1);     Phaa=allP(:,2);     Plna=allP(:,3);     Plca=allP(:,4); 
Paop=allP(:,5);    Prula=allP(:,6);    Prica=allP(:,7);    Plica=allP(:,8);
Plula=allP(:,9);   Psap=allP(:,10);    Prsv=allP(:,11);    Prijv=allP(:,12);
Plijv=allP(:,13);  Plsv=allP(:,14);    Psv=allP(:,15);     Pvc=allP(:,16); 
Pra=allP(:,17);    Prv=allP(:,18);     Prpap=allP(:,19);   Plpap=allP(:,20); 
Prpad=allP(:,21);  Plpad=allP(:,22);   Prpv=allP(:,23);    Plpv=allP(:,24);   
Pla=allP(:,25);

t=0:step:Tall;t=t';
if length(t)<length(Plv)
   t=[t; t(end)+step]; 
end
 
%% Changes in total blood volume during the simulation
i=1;
for i = 1:size(allV,1) 
    
    totalV(i) = sum(allV(i,:));
  
end

%% Changes of blood volume in pulmonary circulation
for i = 1:size(allV,1)
 
   totalVP(i) = sum(allV(i,19:24));

end

%% Changes of blood volume in systemic circulation
for i = 1:size(allV,1)
 
   totalVS(i) = sum(allV(i,2:16));

end

%% Changes of blood volume in the four heart chambers
for i = 1:size(allV,1)
 
   totalVH(i) = allV(i,1)+allV(i,17)+allV(i,18)+allV(i,25);

end
