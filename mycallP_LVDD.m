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

function P=mycallP_LVDD(Vx,t,Fcon,cham,BN)

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
    k8=0.0012;
    Ees(2)=Ees(2)+k8*BN;
elseif cham==3 %---LA
    k9=0.0004;
    k10=0.0008;
    k11=0.0000023;
    Ees(3)=Ees(3)+k9*BN;
    M0(3)=M0(3)+k10*BN;
    lamda(3)=lamda(3)+k11*BN;
end

% The values of parameters in  activation function of four heart chambers
if cham==3  %LA

    ai_0=[ 0      0     0.9    0      0      0      0       0      0     0      ];
    bi=[ 0.12   0.09   0.038  0.07   0.09   0.05    0.04   0.08   0.1    0.04   ];
    ci=[0.005   0.08   0.14   0.25   0.31   0.375   0.45   0.62   0.784  0.7845 ];
    k_a1=0.00044;
    k_a2=0.000188;
    k_a3=-0.000575;
    k_a4=0.0001;
    k_a5=0.0003;
    k_a6=0.000278;
    k_a7=0.000491;
    k_a8=0.000265;
    k_a9=0.000492;
    k_a10=0.000045; 
    ai(1)=ai_0(1,1)+k_a1*BN;
    ai(2)=ai_0(1,2)+k_a2*BN;
    ai(3)=ai_0(1,3)+k_a3*BN;
    ai(4)=ai_0(1,4)+k_a4*BN;
    ai(5)=ai_0(1,5)+k_a5*BN;
    ai(6)=ai_0(1,6)+k_a6*BN;
    ai(7)=ai_0(1,7)+k_a7*BN;
    ai(8)=ai_0(1,8)+k_a8*BN;
    ai(9)=ai_0(1,9)+k_a9*BN;
    ai(10)=ai_0(1,10)+k_a10*BN;
    
else    %LV,RV,RA
   %||---------------ventricle-------------||-----RA---||    
    Ai=[0.3      0.35      0.5        0.55        0.9   ];
    Bi=[0.045    0.035     0.037      0.036       0.038 ];
    Ci=[0.175    0.23      0.275      0.3         0.025 ]+0.1;
end

%% Modify the relative time parameters of atriums and ventricles
if cham==1  %---LV
    Ci(1:4)=Ci(1:4);
elseif cham==2 %---RV
    Ci(1:4)=Ci(1:4)+0.013;
elseif cham==4 %---RA
    Ci(5)=Ci(5);
end

  %% The paratemers of Neuromodulation
  amin=-2;bmin=0.7;Ka=7;Kb=0.5;
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
    for i=1:10
        En=En+ai(i)*exp(-0.5*(((t-ci(i))/bi(i)).^2));
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