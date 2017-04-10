
function dy = NovTyson95eqnse(t,y,k)
% sytem of differential equations based on Novak & Tyon, J. Theor. Biol.,
% 1995

global k1AA;
%global mu;
global kdeg;
global ksynscale;
global kt;
global V;

% synthesis/growth rate constants - kt, kp, ksyn, mu

evarname;

%dy = zeros(18,1); 
dy = zeros(17,1);

% ks - k1AA, k3, kINH, kCAK, ka, kbs, kf1, ke1, kf2,
% ke2, kt, kx, kh, kg, kp, ko, ki, kj, kc, kd, kk, kl, ksyn, kdeg
%k1AA = 0.003;   
%k3 = k3scale*10;     kINH = 0.01;       kCAK = 1;   ka = 0.5;     kbs = 0.2;
%kf1 = 1;        ke1 = 1;     kf2 = 1;           ke2 = 1;    
%kx = 0.1;       kh = 2;      
%kg = 0.2;          kp = 2;     ko = 0.2;
%ki = 0.2;       kj = 0.2;    kc = 0.1;          kd = 0.05;  kk = 1;
%kl = 0.5;       ksyn = ksynscale*0.04; 

%k1AA = k(1);

k3 = k(2);  kINH = k(3);    kCAK = k(4);    ka = k(5);
kbs = k(6);     kf1 = k(7); ke1 = k(8);     kf2 = k(9);     ke2 = k(10);
kx = k(11);     %kh = k(12); 
%kg = k(13);     
kp = ksynscale*k(12);     ko = k(13);
ki = k(14);     kj = k(15); kc = k(16);     kd = k(17);     kk = k(18);
kl = k(19);     ksyn = ksynscale*k(20);

% Ks - Km, Ka (KMA), Kb (KMB), Kf2, Ke2, Kh, Kg, Ki, Kj, Kc, Kd, Kk, Kl
%Km = 1;     Ka = 0.1;     Kb = 0.1;     Kf2 = 0.15;    Ke2 = 0.15;    Kh = 0.05;
%Kg = 0.05;  Ki = 0.01;    Kj = 0.01;    Kc = 0.01;     Kd = 0.01;     Kk = 0.02;
%Kl = 0.02;

Km = k(21);     Ka = k(22); Kb = k(23); Kf2 = k(24);    Ke2 = k(25);
%Kh = k(28);
%Kg = k(29);
Ki = k(26); Kj = k(27);     Kc = k(28);
Kd = k(29);     Kk = k(30); Kl = k(31);

% Vs - V25', V25'', Vwee', Vwee'', V2', V2'', Vmik', Vmik''
%V25pr = 0.045;  V25dpr = 0.45;     Vweepr = 0.0333;         Vweedpr = 0.75;
%V2pr = 0.025;   V2dpr = 0.5;       Vmikpr = 0.02;           Vmikdpr = 0.2;

%F = 15;    Fb = 2;

V25pr = k(32);  V25dpr = k(33);     Vweepr = k(34);     Vweedpr = k(35);
V2pr = k(36);   V2dpr = k(37);      Vmikpr = k(38);     Vmikdpr = k(39);

F = k(40);  Fb = k(41);     mu = ksynscale*0.00385;

% HOW MANY VARIABLES? SOME OF THESE ARE CONSTANTS/HARD-CODED
% y(1) - Cdc13free
% y(2) - MPFphos0
% y(3) - MPFinactive
% y(4) - MPFphos2
% y(5) - MPFactive
                    % y(6) - Cdc2free
% y(6) - Cdc25active
% y(7) - Wee1active
% y(8) - Wee1phosr
% y(9) - Wee1phosl
% y(10) - X
% y(11) - W
% y(12) - Mik1
% y(13) - IEphos
% y(14) - UbE
% y(15) - TAphos
% y(16) - Cdc25total
% y(17) - DNA
% y(18) - mass

PK = 0.965;         Nim1 = 0.258;
UbEtotal = 1;       Wee1total = 1;      Mik1total = 1;
%Wtotal = 1;
IEtotal = 1;        Tatotal = 1;
Cdc2total = 1;

% DERIVED QUANTITIES
Wee1phos2 = Wee1total - y(nWee1active) - y(nWee1phosl) - y(nWee1phosr);
k25 = V25pr*(y(nCdc25total) - y(nCdc25active)) + V25dpr*y(nCdc25active);
kwee = Vweepr*(Wee1total - y(nWee1active)) + Vweedpr*y(nWee1active);
% ADDED SCALAR MULTIPLE OF 0.85 TO CDC13 DEGRADATION RATE FOR OPTIMIZATION
% - 020717
%k2 = 0.85*(V2pr*(UbEtotal - y(nUbE)) + V2dpr*y(nUbE));
k2 = (V2pr*(UbEtotal - y(nUbE)) + V2dpr*y(nUbE));
kmik = Vmikpr*(Mik1total - y(nMik1)) + Vmikdpr*y(nMik1);

% THERE ARE 19 RATE EQUATIONS
% dy(1) - change in free Cdc13
% dy(2) - change in MPFphos0
% dy(3) - change in MPFinactive
% dy(4) - change in MPFphos2
% dy(5) - change in MPFactive
                                    % dy(6) - change in free Cdc2
% dy(6) - change in active Cdc25
% dy(7) - change in active Wee1
% dy(8) - change in "right phosphorylated" Wee1
% dy(9) - change in "left phosphorylated" Wee1
                                    % dy(10) - change in X
% dy(10) - change in W
% dy(11) - change in Mik1
% dy(12) - change in IEphos
% dy(13) - change in UbE
% dy(14) - change in TAphos
% dy(15) - change in total Cdc25
% dy(16) - change in DNA content
% dy(17) - change in mass
% free cyclin, Cdc13
%dy(1) = k1AA - k2*y(nCdc13free) - k3*y(nCdc13free)*y(nCdc2free);

% REMOVED SIZE DEPENDENCE OF CDC13 SYNTHESIS
dy(1) = k1AA - k2*y(nCdc13free) - k3*y(nCdc13free)*(Cdc2total-y(nMPFphos0)-y(nMPFinactive)-y(nMPFphos2)-y(nMPFactive));
%dy(2) = kINH*y(nMPFactive) - ((kwee/(Km+y(nMPFphos0))) + kmik + kCAK + k2)*y(nMPFphos0) + k25*y(nMPFinactive) + k3*y(nCdc13free)*y(nCdc2free);
dy(2) = kINH*y(nMPFactive) - ((kwee/(Km+y(nMPFphos0))) + kmik + kCAK + k2)*y(nMPFphos0) + k25*y(nMPFinactive) + k3*y(nCdc13free)*(Cdc2total-y(nMPFphos0)-y(nMPFinactive)-y(nMPFphos2)-y(nMPFactive));

dy(3) = (kwee*y(nMPFphos0)/(Km+y(nMPFphos0))) + kmik*y(nMPFphos0) - (k25 + kCAK + k2)*y(nMPFinactive) + kINH*y(nMPFphos2);

dy(4) = (kwee*y(nMPFactive)/(Km+y(nMPFactive))) + kmik*y(nMPFactive) - (kINH + k25 + k2)*y(nMPFphos2) + kCAK*y(nMPFinactive);

dy(5) = kCAK*y(nMPFphos0) - (kINH + (kwee/(Km + y(nMPFactive))) + kmik + k2)*y(nMPFactive) + k25*y(nMPFphos2);

%dy(6) = k2*(y(nMPFphos0) + y(nMPFinactive) + y(nMPFphos2) + y(nMPFactive)) - k3*y(nCdc13free)*y(nCdc2free);

%dy(6) = ((ka*y(nMPFactive)*(y(nCdc25total)-y(nCdc25active)))/(Ka + y(nCdc25total) - y(nCdc25active))) - ((kbs*y(nCdc25active))/(Kb + y(nCdc25active)));
dy(6) = ((ka*y(nMPFactive)*(y(nCdc25total)-y(nCdc25active)))/(Ka + y(nCdc25total) - y(nCdc25active))) - ((kbs + (Fb*y(nW)))*y(nCdc25active))/(Kb + y(nCdc25active));
dy(7) = kf1*y(nWee1phosr) - ke1*y(nMPFactive)*y(nWee1active) + ((kf2*y(nWee1phosl))/(Kf2 + y(nWee1phosl))) - ((ke2*y(nWee1active))/(Ke2 + y(nWee1active)))*(PK+Nim1);
dy(8) = ke1*y(nMPFactive)*y(nWee1active) - kf1*y(nWee1phosr) + ((kf2*Wee1phos2)/(Kf2 + Wee1phos2)) - ((ke2*y(nWee1phosr))/(Ke2 + y(nWee1phosr)))*(PK+Nim1);
dy(9) = kf1*Wee1phos2 - ke1*y(nMPFactive)*y(nWee1phosl) - ((kf2*y(nWee1phosl))/(Kf2 + y(nWee1phosl))) + ((ke2*y(nWee1active))/(Ke2 + y(nWee1active)))*(PK+Nim1);
%dy(10) = kt - kx*y(nX);
%dy(10) = ((kh*(Wtotal-y(nW)))/(Kh + Wtotal - y(nW))) - ((kg*y(nW))/(Kg + y(nW)));
dy(10) = kt - kx*y(nW);
dy(11) = kp*y(nW)*(Mik1total - y(nMik1)) - ko*y(nMik1);
dy(12) = ((ki*y(nMPFactive)*(IEtotal - y(nIEphos)))/(Ki + IEtotal - y(nIEphos))) - ((kj*(Tatotal - y(nTaphos))*y(nIEphos))/(Kj + y(nIEphos)));
dy(13) = ((kc*y(nIEphos)*(UbEtotal-y(nUbE)))/(Kc + UbEtotal - y(nUbE))) - ((kd*y(nUbE))/(Kd + y(nUbE)));
dy(14) = ((kk*y(nMPFactive)*(Tatotal - y(nTaphos)))/(Kk + Tatotal - y(nTaphos))) - ((kl*y(nTaphos))/(Kl + y(nTaphos)));
dy(15) = ksyn*y(nmass) - kdeg*y(nCdc25total);
dy(16) = V/(1+(F*y(nTaphos)));
dy(17) = mu*y(nmass);