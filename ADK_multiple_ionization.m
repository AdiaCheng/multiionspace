

% Ionization rates

clear all;

dt = 3.e-15;
Z = 1;
Ui = 13.6;
UH = 13.6;
wa = 4.13e16;
Ea = 510.;
%nstar = 3.69*Z/sqrt(Ui);
nstar = Z*sqrt(UH/Ui);
e = exp(1.);

N = 500;
dE = 1000./N;

for i=1:N

El(i) = 0.1 + (i-1)*dE;
    
%Bruhwiler, PoP, 10, 2022 (2003)
% H-like atoms
W_H_ADK(i) = 1.52e15 * 4^nstar*Ui/(nstar*gamma(2*nstar)) * (20.5*Ui^1.5/El(i))^(2*nstar-1) * exp(-6.83*Ui^1.5/El(i)) * dt;

% Complex atoms, H, l = m = 0

W_00_H(i) = wa*sqrt(3)*(e/pi)^1.5 * (Z^2/nstar^4.5) * (4*e*Ea*Z^3/El(i)/nstar^4)^(2*nstar-1.5) * exp(-2*Ea/(3*El(i)) * (Z/nstar)^3)*dt;
end


% Complex atoms, Ar, l = m = 0
U(1) = 15.75962;	
U(2) = 27.62967;	
U(3) = 40.74;	
U(4) = 59.81;	
U(5) = 75.02;	
U(6) = 91.009;	
U(7) = 124.323;	
U(8) = 143.460;	
U(9) = 422.45;	
U(10) = 478.69;	
U(11) = 538.96;	
U(12) = 618.26;
U(13) = 686.10;	
U(14) = 755.74;	
U(15) = 854.77;	
U(16) = 918.03;	
U(17) = 4120.8857;	
U(18) = 4426.2296;

Z=18;
for i=1:N
El(i) = 0.1 + (i-1)*dE;
%%% U(1) = 15.75962 eV
Z=1;
nstar = Z*sqrt(UH/U(1));
l = 0; m = 0;
W_00_Ar1(i) = wa*sqrt(3)*(e/pi)^1.5 * (Z^2/nstar^4.5) * (4*e*Ea*Z^3/El(i)/nstar^4)^(2*nstar-1.5) * exp(-2*Ea/(3*El(i)) * (Z/nstar)^3)*dt;
l = 1; m = 0;
C_10_Ar1(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
W_10_Ar1(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_10_Ar1(i) * dt;
l = 1; m = 1;
C_11_Ar1(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
W_11_Ar1(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_11_Ar1(i) * dt;
% l = 2; m = 0;
% C_20_Ar1(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
% W_20_Ar1(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_20_Ar1(i) * dt;
% l = 2; m = 1;
% C_21_Ar1(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
% W_21_Ar1(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_21_Ar1(i) * dt;
% l = 1; m = 1;
% C_22_Ar1(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
% W_22_Ar1(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_22_Ar1(i) * dt;


%%% U(3) = 40.74 eV
Z=3;
nstar = Z*sqrt(UH/U(2));
l = 0; m = 0;
W_00_Ar3(i) = wa*sqrt(3)*(e/pi)^1.5 * (Z^2/nstar^4.5) * (4*e*Ea*Z^3/El(i)/nstar^4)^(2*nstar-1.5) * exp(-2*Ea/(3*El(i)) * (Z/nstar)^3)*dt;
l = 1; m = 0;
C_10_Ar3(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
W_10_Ar3(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_10_Ar3(i) * dt;
l = 1; m = 1;
C_11_Ar3(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
W_11_Ar3(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_11_Ar3(i) * dt;
% l = 2; m = 0;
% C_20_Ar3(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
% W_20_Ar3(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_20_Ar3(i) * dt;
% l = 2; m = 1;
% C_21_Ar3(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
% W_21_Ar3(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_21_Ar3(i) * dt;
% l = 2; m = 2;
% C_22_Ar3(i) = (2*l+1)*factorial(l+m) / (2*pi*nstar*2^m * factorial(m)*factorial(l-m));
% W_22_Ar3(i) = wa*sqrt(3*nstar^3/pi/Z^3 *El(i)/Ea) * Z^2/(2*nstar^2)* (2*e/nstar)^(2*nstar) * exp(-2/3*Ea/El(i)*(Z/nstar)^3) * (2*Ea/El(i)*Z^3/nstar^3)^(2*nstar-m-1) * C_22_Ar3(i) * dt;

end


%plot(El,W_H_ADK, El,W_00_H,  El,W_00_Ar1,  El,W_00_Ar3, 'linewidth',3)
semilogy(El,W_H_ADK,  El,W_00_H,  El,W_00_Ar1,  El,W_00_Ar3, 'linewidth',3)
%loglog(El,W_H_ADK,  El,W_00_H,  El,W_00_Ar1,  El,W_00_Ar3, 'linewidth',3)
set(gca,'fontsize',30);
grid on
xlabel('Laser field, GV/m','FontSize',30);
ylabel('Ionization Rate','FontSize',30);
legend('ADK H','00 H','00 Ar1','00 Ar3');

figure
semilogy(El,W_00_Ar1,  El,W_10_Ar1, El,W_11_Ar1, 'linewidth',3)
set(gca,'fontsize',30);
grid on
xlabel('Laser field, GV/m','FontSize',30);
ylabel('Ionization Rate','FontSize',30);
legend('00 Ar1','10 Ar1', '11 Ar1');

figure
semilogy(El,W_00_Ar3,  El,W_10_Ar3, El,W_11_Ar3, 'linewidth',3)
set(gca,'fontsize',30);
grid on
xlabel('Laser field, GV/m','FontSize',30);
ylabel('Ionization Rate','FontSize',30);
legend('00 Ar3','10 Ar3', '11 Ar3');
% 
% k1 = exp(-2/3*Ea./El*(Z/nstar)^3);
% k2 = sqrt(3*nstar^3/pi/Z^3 *El/Ea);
% k3 = (2*Ea./El*Z^3/nstar^3).^(2*nstar-m-1);
% figure
% semilogy(El,k1, El,k2, El,k3, 'linewidth',3)
% set(gca,'fontsize',30);
% grid on
% xlabel('Laser field, GV/m','FontSize',30);
% %ylabel('Ionization Rate','FontSize',30);
% legend('k1','k2','k3');
% 
% Z=1; l=0; m=0;
% nstar = 3.69*Z/sqrt(Ui);
% k1 = exp(-2/3*Ea./El*(Z/nstar)^3);
% k2 = sqrt(3*nstar^3/pi/Z^3 *El/Ea);
% k3 = (2*Ea./El*Z^3/nstar^3).^(2*nstar-m-1);
% figure
% semilogy(El,k1, El,k2, El,k3, 'linewidth',3)
% set(gca,'fontsize',30);
% grid on
% xlabel('Laser field, GV/m','FontSize',30);
% %ylabel('Ionization Rate','FontSize',30);
% legend('k1','k2','k3');
% title('Hydrogen, Z=1, l=m=0, U_H = 13.6eV');
