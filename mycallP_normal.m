% function P=mycallP(Vx,t,Fcon,cham)
% This function calculates the blood pressure of four heart chambers

%% Inputs:
% Vx is volumle of four heart chambers
% t is time 
% Fcon is sympathetic efferent discharge frequency

%% Outputs: 
% The output is blood pressure
% cham==1, The output is blood pressure of left ventricle
% cham==2, The output is blood pressure of right ventricle
% cham==3, The output is blood pressure of left atrium
% cham==4, The output is blood pressure of right atrium

function P=mycallP(Vx,t,Fcon,cham)

%The values of parameters in P-V relationship of four heart chambers
%----LV         RV          LA          RA   
Ees=[4.3        0.8         0.3        0.3   ];
M0=[1.7         0.67        0.5        0.5   ];
V0=[25          25          20         20    ];
Vd=[40          40          20         20    ];
lamda=[0.015   0.015     0.025       0.025   ];

% The values of parameters in  activation function of four heart chambers
%||---------------ventricle------------||----atrium---||
Ai=[0.3      0.35      0.5        0.55        0.9   ];
Bi=[0.045    0.035     0.037      0.036       0.038 ];
Ci=[0.175    0.23      0.275      0.3         0.025 ]+0.1;

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
  amin=-2;bmin=0.7;Ka=7;Kb=0.5;
  a=amin+Ka*Fcon;
  b=bmin+Kb*Fcon;

if cham==1 %  left ventricle
    En=0;
    for i=1:4
        En=En+Ai(i)*exp(-0.5*(((b.*t-Ci(i))/Bi(i)).^2));
    end
    Pes=a.*Ees(cham).*(Vx-Vd(cham));
    Ped=M0(cham)*abs(exp(    lamda(cham)  *   (   Vx  -  V0(cham)   )      )-1);
    P=En.*Pes+(1-En).*Ped;
    
else if cham==2   % right ventricle      
    En=0;
    for i=1:4
        En=En+Ai(i)*exp(-0.5*(((b.*t-Ci(i))/Bi(i)).^2));
    end
    Pes=a.*Ees(cham).*(Vx-Vd(cham));        
    Ped=M0(cham)*abs(exp(    lamda(cham) *   (   Vx  -  V0(cham)   )      )-1);
    P=En.*Pes+(1-En).*Ped;
    else      % left and right atriums
    En=0;
    for i=5:5
        En=En+Ai(i)*exp(-0.5*(((t-Ci(i))/Bi(i)).^2));
    end
    Pes=Ees(cham)*(Vx-Vd(cham));
    Ped=M0(cham)*abs(exp(lamda(cham)*(Vx-V0(cham)))-1);
    P=En.*Pes+(1-En).*Ped;
    end
end