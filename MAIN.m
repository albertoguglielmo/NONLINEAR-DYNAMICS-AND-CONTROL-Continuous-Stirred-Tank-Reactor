close all
clear

%per far paritre lo script si possono decidere:
u_start_lin=0;
steps=1;
u_end_lin=40;
load('Q4.mat') %or load('Q5.mat')
load('bin_u.mat')
load('bin_x1.mat')
load('bin_x2.mat')


%% parameters

% Reactor Parameters
V = 50; % Reactor volume (l)
F= 50; % Inlet volumetric flow rate to the reactor (l/min)
C_A0 = 1; % Feed concentration of component A (mole/l)
K_0 = 7.8e10; % Pre-exponential factor (l/min)
E_over_R = 8567; % Activation energy in the Arrhenius equation (Cal/mole)
R = 1.987; % Universal gas constant (Cal/mole.K), assuming typical value
rho = 900; % Density of the inlet and outlet stream (g/l)
C_p = 0.329; % Heat capacity of inlet and outlet stream (Cal/g.K)
T_in = 350; % Inlet stream temperature (K)
H_r = -5e4; % Heat of reaction (Cal/mole)
UA = 5e4; % Heat transfer term (Cal/min.K)
rho_c = 1000; % Density of the coolant water in the jacket (g/l)
T_amb=293; %Temperatura ambiente (K)
xi_0 = [0; 0];  
Tabella = [];    
n = 1;   

for in = u_start_lin:steps:u_end_lin

    xi_0 = [0; 0];
    handle = @(x) [(F/V) *(C_A0 - x(1)) - K_0 * x(1) * exp(-E_over_R / x(2));...
   (F/V) * (T_in - x(2))- (K_0 *x(1)*  H_r) / (rho * C_p) * exp(-E_over_R / x(2)) - (UA / (rho * V * C_p)) * ( x(2)- (T_amb)-in)];

    xi_eq = fsolve(handle, xi_0);
    
    Tabella(n, :) = [in, xi_eq(1), xi_eq(2)];
    
    n = n + 1;
end




V = 50; % Reactor volume (l)
F= 50; % Inlet volumetric flow rate to the reactor (l/min)
C_A0 = 1; % Feed concentration of component A (mole/l)
K_0 = 7.8e10; % Pre-exponential factor (l/min)
E_over_R = 8567; % Activation energy in the Arrhenius equation (Cal/mole)
R = 1.987; % Universal gas constant (Cal/mole.K), assuming typical value
rho = 900; % Density of the inlet and outlet stream (g/l)
C_p = 0.329; % Heat capacity of inlet and outlet stream (Cal/g.K)
T_in = 350; % Inlet stream temperature (K)
H_r = -5e4; % Heat of reaction (Cal/mole)
UA = 5e4; % Heat transfer term (Cal/min.K)
T_amb=293; %Temperatura ambiente (K)

for in=1:length(Tabella)

u=Tabella(in,1);
C_A =   Tabella(in,2)   ;         % (mole/l)
T =    Tabella(in,3)   ;         % (K)

A=[-(F/V)-K_0*exp(-E_over_R / T)  ,- K_0 * C_A * exp(-E_over_R / T)/T^2*(-E_over_R)
       - (K_0 *  H_r) / (rho * C_p) * exp(-E_over_R / T) ,-(F/V) - (K_0 * C_A *  H_r) / (rho * C_p) * exp(-E_over_R / T)/T^2*(-E_over_R)- (UA / (rho * V * C_p))];
B = [0; 3.377];
C=[1 0] ;

zita=0.7;
wn=2.3;


K2=(A(1,1)+A(2,2)+2*zita*wn)/B(2);
K1=(+wn^2+A(1,2)*A(2,1)-A(1,1)*A(2,2)+K2*B(2)*A(1,1))/(B(2)*A(1,2));
Kr=-1/(C*inv(A-B*[K1 K2])*B);




Tabella(in,4)=K1;
Tabella(in,5)=K2;
Tabella(in,6)=Kr;
end








%per la linearizzazione
C_A =   0.554   ;         % Concentration of component A in outlet stream (mole/l)
T =    3.386561081008934e+02   ;         % Temperature of the reactants in the reactor (K)
A=[-(F/V)-K_0*exp(-E_over_R / T)  ,- K_0 * C_A * exp(-E_over_R / T)/T^2*(-E_over_R)
       - (K_0 *  H_r) / (rho * C_p) * exp(-E_over_R / T) ,-(F/V) - (K_0 * C_A *  H_r) / (rho * C_p) * exp(-E_over_R / T)/T^2*(-E_over_R)- (UA / (rho * V * C_p))];
B = [0; 3.377];
C=[1 0] ;

K2=(A(1,1)+A(2,2)+2*zita*wn)/B(2);
K1=(+wn^2+A(1,2)*A(2,1)-A(1,1)*A(2,2)+K2*B(2)*A(1,1))/(B(2)*A(1,2));
Kr=-1/(C*inv(A-B*[K1 K2])*B);




%% functions
function xi_dot = react(x)

% Reactor Parameters
V = 50; % Reactor volume (l)
F= 50; % Inlet volumetric flow rate to the reactor (l/min)
C_A0 = 1; % Feed concentration of component A (mole/l)
K_0 = 7.8e10; % Pre-exponential factor (l/min)
E_over_R = 8567; % Activation energy in the Arrhenius equation (Cal/mole)
R = 1.987; % Universal gas constant (Cal/mole.K), assuming typical value
rho = 900; % Density of the inlet and outlet stream (g/l)
C_p = 0.329; % Heat capacity of inlet and outlet stream (Cal/g.K)
T_in = 350; % Inlet stream temperature (K)
H_r = -5e4; % Heat of reaction (Cal/mole)
UA = 5e4; % Heat transfer term (Cal/min.K)
rho_c = 1000; % Density of the coolant water in the jacket (g/l)
T_amb=293; %Temperatura ambiente (K)

xi_dot = [(F/V) *(C_A0 - x(1)) - K_0 * x(1) * exp(-E_over_R / x(2));...
   (F/V) * (T_in - x(2))- (K_0 *x(1)*  H_r) / (rho * C_p) * exp(-E_over_R / x(2)) - (UA / (rho * V * C_p)) * ( x(2)- (T_amb)-20)];


end
