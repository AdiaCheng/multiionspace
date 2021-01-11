
% source parameters
ep0 = 8.8542e-12; % (F/m) vacuume permitivoty
mu0 = pi*4e-7; % (H/m) vacuume permeability
c = 2.998e8;        % (m/s)
l0 = 10.6e-6; %(m) wavelength
f0 = c/l0; %(Hz) frequency
tau0 = 1/f0; %(s) period

dz = l0/20; %(m) grid resolution
n = 0;  % time index
dt = tau0/50; % (s) time resolution
tc = 25*tau0/2; % (s) center at # cycles of signal
sigma = tc/(3*sqrt(2*log(2)));
% gt = exp(-(n*dt-tc)^2/(2*sigma^2)); % gaussian pulse
% mt = sin(2*pi*f0*n*dt); % modulation

% field initialization
Distance = 3e-4; %(m)
E0 = 500;   % (GV)
M = ceil(Distance/dz);
z = linspace(0,Distance, M); % (m) z coordinate
Hy = zeros(1,M);
Ex = zeros(1,M+1);
nppc = 40;  % number of particles per cell.
n0 = 1.7e16; % (g/m^3) initail neutral particle density (nc from Chen.-JCP)
energyAr = [15.75962, 27.62967, 40.74, 59.81, 75.02, 91.009, 124.323, 143.460, 422.45, 478.69, 538.96, 618.26, 686.10, 755.74, 854.77, 918.03, 4120.8857, 4426.2296]; % (eV)

% result initialization
Ar = zeros(length(energyAr)+1,M+1);
newAr = zeros(length(energyAr)+1, M+1);
Ar(1,:) = ones(1,M+1) * nppc; % one row of neutral
P = zeros(length(energyAr),M+1);
coefn = zeros(length(energyAr),1);
coef1 = zeros(length(energyAr),1);
coef2 = zeros(length(energyAr),M+1);
W = zeros(length(energyAr),M+1);



% insert source and update field
for n = 1:1300
    
    ndt = n*dt;
    mt = sin(2*pi*f0*ndt); % modulation
    gt = exp(-(ndt-tc)^2/(2*sigma^2)); % gaussian pulse
    source = mt*gt;
    Ex(1) = source;
    Hy = Hy-dt./mu0.*diff(Ex)/dz;
    Ex(2:M) = Ex(2:M)-dt./ep0.*diff(Hy/dz);
    
    Ea = abs(E0*Ex);
    % independent from Ea; same along entire z
    for ion = 1:length(energyAr)
        energy = energyAr(ion);
        coefn(ion) = 3.69*18/sqrt(energy);
        coef1(ion) = 1.52e+15*(4^coefn(ion))*energy/(coefn(ion)*gamma(2*coefn(ion)));
    end
    
    % only calculate when Ea>0
    update = 0;
    for grid = 1:M+1
        if Ea(grid)-eps>0
            update = grid; % records the maximum grid number reached this condition

            for ion = 1:length(energyAr)
                energy = energyAr(ion);
                coef2(ion,grid) = 20.5* (energy^1.5) /Ea(grid);
                W(ion,grid) = coef1(ion)* coef2(ion,grid)*exp(-6.83*energy^1.5/Ea(grid));
            end
            
            
            %%%%% ionization solution
%             % --- based on explicit solution of differential equation ---
%  
%             % n^0            
%             P(1,grid) = (1-exp(-W(1,grid)*dt));
%             newAr(1,grid) = Ar(1,grid)*(1-P(1,grid));
%             
%             % n^1+ ~ n^17+
%             for ion = 2:length(energyAr)
%                 % n^i+ corresponds to energyAr(i+1)
%                 if W(ion,grid) < eps
%                     newAr(ion,grid) = W(ion-1,grid)*Ar(ion-1,grid)*dt;
%                 else
%                     P(ion,grid) = (1-exp(-W(ion,grid)*dt));
%                     %add = W(ion-1,grid)/W(ion,grid) * Ar(ion-1,grid);
%                     %newAr(ion,grid) = (1-P(ion,grid))*Ar(ion,grid) + P(ion,grid)*add;
%                     newAr(ion,grid) =(W(ion,grid)*Ar(ion,grid)-W(ion-1,grid)*Ar(ion-1,grid))/W(ion,grid) * (1-P(ion,grid)) + W(ion-1,grid)/W(ion,grid) * Ar(ion-1,grid);
%                 end                
%             end
%             
%             % n^18+ (in shorter test version, last ion)
%             ion = ion+1;
%             newAr(ion,grid) = Ar(ion,grid) + W(ion-1,grid)*Ar(ion-1,grid)*dt;
           
            
            
            % --- based on 4th order RK approximation of ODE ---
            
            % n^0
            f{1} = @(t, y) -W(1,grid)*y;
            newAr(1,grid) = RK4(f{1}, ndt, dt, Ar(1,grid));
            
            % n^1+ ~ n^17+
            for ion = 2:length(energyAr)
                f{ion} = @(t, y) W(ion-1,grid)*Ar(ion-1,grid)-W(ion,grid)*y;
                newAr(ion,grid) = RK4(f{ion}, ndt, dt, Ar(ion,grid));
            end
            
            % n^18+
            f{19} = @(t,y) W(18,grid)*Ar(18,grid); % Ar(18)~ n^17+
            newAr(19,grid) = RK4(f{19}, ndt, dt, Ar(19,grid));
            
            
        end
    end
    
    % %%% check propagation
    %    plot(Ex)
    %    legend(int2str(n))
    %    pause(1/1000)
    %    surf(coef2)
    %   surf(W)
    %   zlabel('W (s^{-1})')
    %   pause(1/1000)
    
    
     Ar(:,1:update) = newAr(:,1:update);
     %%%% check ionization
     surf(Ar)
     zlabel('no. macro particles')
     pause(1/1000)



    
    
    
end

xlabel('longitudinal cell indices')
ylabel('ionization level')


figure
subplot(5,1,1)
plot(Ex)
ylabel('normalized E field')
subplot(5,1,2)
plot(Ar(2,:))
ylabel('Ar^+')
subplot(5,1,3)
plot(Ar(3,:))
ylabel('Ar^{2+}')
subplot(5,1,4)
plot(Ar(4,:))
ylabel('Ar^{3+}')
subplot(5,1,5)
plot(Ar(12,:))
ylabel('Ar^{11+}')
xlabel('longitudinal cell indices')
title('no. macroparticles')



function y = RK4(f, t, h, y)
    k1 = f(t,y);
    k2 = f(t+h/2, y+h/2*k1);
    k3 = f(t+h/2, y+h/2*k2);
    k4 = f(t+h, y+h*k3);
    y = y + h/6 * (k1+2*k2+2*k3+k4);
end




