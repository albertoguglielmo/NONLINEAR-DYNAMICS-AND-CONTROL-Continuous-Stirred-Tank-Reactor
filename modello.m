function [A,B,C,D] = modello (C_A,T)

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


    % Calcolo delle derivate
    %x1dot = (F/V) *(C_A0 - C_A) - K_0 * C_A * exp(-E_over_R / T); %Cadot

A=[-(F/V)-K_0*exp(-E_over_R / T)  ,-(-E_over_R)*K_0 * C_A * exp(-E_over_R / T)/T^2
       - (K_0 *  H_r) / (rho * C_p) * exp(-E_over_R / T) ,-(F/V) - (K_0 * C_A *  H_r) / (rho * C_p) * exp(-E_over_R / T)*(-E_over_R)/T^2- (UA / (rho * V * C_p))];
B = [0; 3.377];
C=[1 0] ;
D=0;

end