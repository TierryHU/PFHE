function f = funcNs2019(x,D)
% Plate fin heatexchenger , Ns NTU
% Zhang, L., Hu, T., Yang, Z. et al. 
%Elite and dynamic opposite learning enhanced sine cosine algorithm for application to plat-fin heat exchangers design problem. 
%Neural Comput & Applic (2021). https://doi.org/10.1007/s00521-021-05963-2
% Coded by Tianyu Hu (Thierry)
% Modified on 21-02-2020
% Greek letters execute LaTeX language
%h-hot,c-cold

Lh = x(1);%[0.1,1]             
Lc = x(2);%1000kw
H = x(3);%m,mm[0.002,0.01]
t = x(4);%m,mm[0.0001,0.0002]
n = x(5);%[100,1000]
l = x(6);%m,mm[0.001,0.01]
Nh = x(7); 

% mh   = 0.8962;   % Mass flow rate, m (kg s^-1)
% mc   = 0.8296; 
% Th  = 513;      % Inlet temperature, T1 (K)
% Tc  = 277; 
% Ph  = 10^5;     % Inlet pressure, P1 (Pa)
% Pc  = 10^5;
% Cph  = 1017.7;   % Specific heat, Cp (J kg^-1 K^-1)
% Cpc  = 1011.8;
% rhoh = 0.8196;   % Density, q (kg m^-3)
% rhoc = 0.9385;
% 
% muh = 0.0000241;
% muc = 0.00002182;
% Prh  = 0.6878;    % Prandtl number, Pr
% Prc  = 0.6954;
% Ch = mh*Cph;
% Cmax = Ch;
% Cc = mc*Cpc;
% Cmin = Cc;
% 
% Q    = 159.99;       % Heat duty of the exchanger, Q (kW)



%Case 1
mh = 1.66;
mc = 2;
Th = 900+273.15;%K
Tc = 200+273.15;
Ph = 160000;%Pa
Pc = 200000;
Cph = 1122;
Cpc = 1073;
rhoh = 0.6296;%paper
rhoc = 0.9638;
muh = 0.0000401;
muc = 0.0000336;
Prh = 0.731;
Prc = 0.694;
Ch = mh*Cph;
Cc = mc*Cpc;
Cmax = Cc;
Cmin = Ch;
%Q=1070kw
Hh = H;
Hc = H;

nh = n;
nc = n;

s = (1/n-t);%(34)
%Dh = (2*(s-t)*(H-t))/(s+(H-t)+((H-t)*t)/l);%(33)
%dh = (2*(s-t)*(H-t))/((s+(H-t))+((H-t)*t)/(l)); %(33)2009
dh = (4*s*(H-t)*l)/(2*(s*l+((H-t)*l)+(t*(H-t)))+t*s);%2019
Dhh = dh;
Dhc = dh;

Nc = Nh+1;

Cr = Cmin/Cmax;

Affh = (Hh-t)*(1-(nh*t))*Lc*Nh;%(20)
Affc = (Hc-t)*(1-(nc*t))*Lh*Nc;%(21)

Gh = mh/Affh; %(m/Aff)
Gc = mc/Affc;

%St = h/G*Cp;
%Sta = h/Ga*Cpa;
%Stb = h/Gb*Cpb;

%Rea = (Ga*Dha)/mua; %
%Reb = (Gb*Dhb)/mub; %
Reh = (mh*dh)/(Affh*muh); % (32)
Rec = (mc*dh)/(Affc*muc); %

%For laminar flow (Re＜1500)
if Reh>1500

        jh = 0.21*(Reh^(-0.4))*((l/dh)^(-0.24))*(t/dh)^(0.02);%(30)
        fh = 1.12*(Reh^(-0.36))*((l/dh)^(-0.65))*(t/dh)^(0.17);%(31)
else
            %<=1500
        jh = 0.53*(Reh^(-0.5))*((l/dh)^(-0.15))*((s/H)-t)^(-0.14);%(28)
       fh = 8.12*(Reh^(-0.74))*((l/dh)^(-0.41))*((s/H)-t)^(-0.02);%(29)
        
end
    
if Rec>1500
        jc = 0.21*(Rec^(-0.4))*((l/dh)^(-0.24))*(t/dh)^(0.02);%(30)
        fc = 1.12*(Rec^(-0.36))*((l/dh)^(-0.65))*(t/dh)^(0.17);%(31)
else
        jc = 0.21*(Rec^(-0.4))*((l/dh)^(-0.24))*(t/dh)^(0.02);
        fc = 1.12*(Rec^(-0.36))*((l/dh)^(-0.65))*(t/dh)^(0.17);
        
end
%deltaPa = (2*fa*(ma^2)*La)/(rhob*Dha*(Lb^2)*(Na^2)*((Ha-ta)^2)*((1-na*ta)^2));%(26)
%deltaPb = (2*fb*(mb^2)*Lb)/(rhoa*Dhb*(La^2)*(Nb^2)*((Hb-tb)^2)*((1-nb*tb)^2));%(27)
deltaPh = (2*fh*Lh*(Gh^2))/(rhoh*Dhh);%deltaP应小于10^5
deltaPc = (2*fc*Lc*(Gc^2))/(rhoc*Dhc);

Ah = Lh*Lc*Nh*(1+2*nh*(Hh-t));%(22)
Ac = Lh*Lc*Nc*(1+2*nc*(Hc-t));%(23)
%AHT = Aa+Ab;%(24)
%NTU = (ja*Cpa*(Pra^(-2/3))*ma*Aa)/(Cmin*Affa)+(jb*Cpb*(Prb^(-2/3))*mb*Ab)/(Cmin*Affb);%(19)
NTU = (1/Cmin)*((jh*Cph*(Prh^(-(2/3)))*mh*Ah)/Affh+(jc*Cpc*(Prc^(-(2/3)))*mc*Ac)/Affc);
%epsilon = 1-exp((1/Cr)*(NTU^0.22)*(exp(-Cr*(NTU^0.78))-1));%(17)
epsilon = 1-exp((1/Cr)*(NTU^0.22)*(exp((-Cr)*(NTU^(0.78)))-1));

Q = epsilon*Cmin*(Th-Tc);
% Qx = (160000 - Q);
Qx = (1070000 - Q);
Qx = abs(Qx);
%Ra = Pa1/(Ta1*rhoa);
%Rb = Pb1/(Tb1*rhob);
%Ra = Ph*(1/rhoh)/Th;
%Rb = Pc*(1/rhoc)/Tc;
%Ns = (Ch/Cmax)*(exp(1-epsilon*(Cmin/Ch)*(1-(Tc1/Th1)))-(Ra/Cph)*exp(1-(deltaPh/Ph1)))+(Cc/Cmax)*(exp(1+epsilon*(Cmin/Cc)*((Th1/Tc1)-1))-(Rb/Cpc)*exp(1-(deltaPc/Pc1)))+(0.01*Qx);
%(16)
%(25)
Retch = (Ph*(1/rhoh))/Th;
Retcc = (Pc*(1/rhoc))/Tc;
%Ns = ((1-epsilon)*(((Tc-Th)^2)/(Tc*Th))+((Retch/Cph)*(deltaPh/Ph))+((Retcc/Cpc)*(deltaPc/Pc)));%160kw
Rh = Ph*(1/rhoh)/Th;
Rc = Pc*(1/rhoc)/Tc;
%Ns = (Ch/Cmax)*(exp(1-(epsilon*(Cmin/Ch)*(1-(Tc/Th))))-(Rh/Cph)*(exp(1-(deltaPh/Ph))))+(Cc/Cmax)*(exp(1+epsilon*(Cmin/Cc)*((Th/Tc)-1))-(Rc/Cpc)*exp(1-(deltaPc/Pc)))+(0.0001*Qx);
Ns = ((1-epsilon)*(((Tc-Th)^2)/(Tc*Th))+((Retch/Cph)*(deltaPh/Ph))+((Retcc/Cpc)*(deltaPc/Pc)))+(0.0001*Qx);%160kw

f = Ns;

