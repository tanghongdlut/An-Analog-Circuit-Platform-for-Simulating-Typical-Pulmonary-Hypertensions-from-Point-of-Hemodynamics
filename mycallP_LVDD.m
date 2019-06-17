% function P=mycallP(Vx,t,Fcon,cham,BN)
% This function calculates the blood pressure of four heart chambers
%% Inputs:
% Vx is volumle of four heart chambers
% t is time 
% Fcon is sympathetic efferent discharge frequency
% BN is the number of cardiac cycles
%% Outputs: 
% The output is blood pressure
% cham==1, The output is blood pressure of left ventricle
% cham==2, The output is blood pressure of right ventricle
% cham==3, The output is blood pressure of left atrium
% cham==4, The output is blood pressure of right atrium

function P=mycallP(Vx,t,Fcon,cham,BN)

%The  values of parameters in P-V relationship of four heart chambers
%----LV         RV          LA          RA   
Ees=[4.3        0.8         0.3        0.3  ];
M0=[1.7        0.67         0.5        0.5  ];
V0=[25          25          20         20   ];
Vd=[40          40          20         20   ];
lamda=[0.015   0.015     0.025       0.025  ];

if cham==1  %---LV
    k3=0.004;
    k4=0.00001;
    M0(1)=M0(1)+k3*BN;
    lamda(1)=lamda(1)+k4*BN;
elseif cham==2 %---RV
    k9=0.0012;
    Ees(2)=Ees(2)+k9*BN;
elseif cham==3 %---LA
    k10=0.0004;
    k11=0.0008;
    k12=0.0000023;
    Ees(3)=Ees(3)+k10*BN;
    M0(3)=M0(3)+k11*BN;
    lamda(3)=lamda(3)+k12*BN;
end

% The values of parameters in  activation function of four heart chambers
if cham==3  %LA
    k_a1=0.0004488;
    k_a2=0.000168;
    k_a3=0.00039;
    k_a4=0.0002;
    k_a5=0.0002;
    k_a6=0.000168;
    k_a7=0.000291;
    k_a8=0.00028;
    k_a9=0.0005129;

    Ai(1)=k_a1*BN;
    Ai(2)=k_a2*BN;
    Ai(3)=0.9-k_a3*BN;
    Ai(4)=k_a4*BN;
    Ai(5)=k_a5*BN;
    Ai(6)=k_a6*BN;
    Ai(7)=k_a7*BN;
    Ai(8)=k_a8*BN;
    Ai(9)=k_a9*BN;
    Bi=[ 0.1    0.05   0.038   0.07   0.1   0.045   0.04    0.07   0.1];
    Ci=[0.005   0.08   0.136   0.25   0.3   0.35    0.42    0.6    0.78];
else    %LV,RV,RA
    Ai=[0.3      0.35      0.5        0.55        0.9   ];
    Bi=[0.045    0.035     0.037      0.036       0.038 ];
    Ci=[0.175    0.23      0.275      0.3         0.025 ]+0.1;
end


%% Modify the relative time parameters of atriums and ventricles
if cham==1  %---LV
    Ci(1:4)=Ci(1:4);
elseif cham==2 %---RV
    Ci(1:4)=Ci(1:4)+0.013;
elseif cham==3 %---LA
    Ci(5)=Ci(5)+0.02;
elseif cham==4 %---RA
    Ci(5)=Ci(5);
end

%% The paratemers of Neuromodulation
  amin=0;bmin=0.7;Ka=3;Kb=0.5;
  a=amin+Ka*Fcon;
  b=bmin+Kb*Fcon;

if cham==1  %LV 
    En=0;
    for i=1:4
        En=En+Ai(i)*exp(-0.5*(((b.*t-Ci(i))/Bi(i)).^2));
    end
    
    Pes=a.*Ees(cham).*(Vx-Vd(cham));
    Ped=M0(cham)*abs(exp(    lamda(cham)  *   (   Vx  -  V0(cham)   )      )-1);
    P=En.*Pes+(1-En).*Ped;
    
else if cham==2   %RV      
    En=0;
    for i=1:4
        En=En+Ai(i)*exp(-0.5*(((b.*t-Ci(i))/Bi(i)).^2));
    end
    
    Pes=a.*Ees(cham).*(Vx-Vd(cham));        
    Ped=M0(cham)*abs(exp(    lamda(cham) *   (   Vx  -  V0(cham)   )      )-1);
    P=En.*Pes+(1-En).*Ped;
    
else if cham==3  %LA
    En=0;
    for i=1:9
        En=En+Ai(i)*exp(-0.5*(((t-Ci(i))/Bi(i)).^2));
    end
    
    Pes=Ees(cham)*(Vx-Vd(cham));
    Ped=M0(cham)*abs(exp(lamda(cham)*(Vx-V0(cham)))-1);
    P=En.*Pes+(1-En).*Ped;
        
    else %RA
    En=0;
    for i=5:5
        En=En+Ai(i)*exp(-0.5*(((t-Ci(i))/Bi(i)).^2));
    end
    
    Pes=Ees(cham)*(Vx-Vd(cham));
    Ped=M0(cham)*abs(exp(lamda(cham)*(Vx-V0(cham)))-1);
    P=En.*Pes+(1-En).*Ped;
    
    end
   end
end