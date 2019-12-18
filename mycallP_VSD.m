% function P=mycallP_VSD(Vx,t,Fcon,cham,BN)
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

function P=mycallP_VSD(Vx,t,Fcon,cham,BN)

%The values of parameters in P-V relationship of four heart chambers
%----LV         RV          LA         RA   
Ees=[4.3        0.8        0.3       0.3   ];
M0=[1.7         0.67       0.5       0.5   ];
V0=[25          25         20        20    ];
Vd=[40          40         20        20    ];
lamda=[0.015   0.015      0.025      0.025 ];

if cham==2 %---RV
    krv_v=0.0012;
    Ees(2)=Ees(2)+krv_v*BN;
elseif cham==3 %---LA
    kla_v= 0.0003;
    Ees(3)=Ees(3)+kla_v*BN;
end

% The values of parameters in  activation function of four heart chambers
if cham==3 %---LA    
    ai_0=[ 0.9      0      0  ];
    bi=[ 0.038    0.05   0.05 ];
    ci=[ 0.145    0.37   0.42 ];
    k_a2=0.000056;
    k_a3=0.000056;
    ai(1,1)= ai_0(1,1);
    ai(1,2)= ai_0(1,2)+k_a2*BN;
    ai(1,3)= ai_0(1,3)+k_a3*BN; 
else
    %||---------------ventricle-------------||---RA--||
    Ai=[0.3      0.35      0.5        0.55      0.9    ];
    Bi=[0.045    0.035     0.037      0.036     0.038  ];
    Ci=[0.175    0.23      0.275      0.3       0.025  ]+0.1;   
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

if cham==1  % left ventricle
    En=0;
    for i=1:4
        En=En+Ai(i)*exp(-0.5*(((b.*t-Ci(i))/Bi(i)).^2));
    end
    
    Pes=a.*Ees(cham).*(Vx-Vd(cham));
    Ped=M0(cham)*abs(exp(    lamda(cham)  *   (   Vx  -  V0(cham)   )      )-1);
    P=En.*Pes+(1-En).*Ped;
    
else if cham==2 % right ventricle        
    En=0;
    for i=1:4
        En=En+Ai(i)*exp(-0.5*(((b.*t-Ci(i))/Bi(i)).^2));
    end
    
    Pes=a.*Ees(cham).*(Vx-Vd(cham));        
    Ped=M0(cham)*abs(exp(    lamda(cham) *   (   Vx  -  V0(cham)   )      )-1);
    P=En.*Pes+(1-En).*Ped;
    
    else if cham==3 % left atrium
    En=0;
    for i=1:3
        En=En+ai(i)*exp(-0.5*(((t-ci(i))/bi(i)).^2));
    end
    
    Pes=Ees(cham)*(Vx-Vd(cham));
    Ped=M0(cham)*abs(exp(lamda(cham)*(Vx-V0(cham)))-1);
    P=En.*Pes+(1-En).*Ped;
        
    else % right atrium
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