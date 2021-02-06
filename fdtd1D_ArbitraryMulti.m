
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
E0 = 400;   % (GV)
M = ceil(Distance/dz);
z = linspace(0,Distance, M); % (m) z coordinate
Hy = zeros(1,M);
Ex = zeros(1,M+1);
nppc = 40;  % number of particles per cell.
n0 = 1.7e16; % (g/m^3) initail neutral particle density (nc from Chen.-JCP)
energyAr = [15.75962, 27.62967, 40.74, 59.81, 75.02, 91.009, 124.323, 143.460, 422.45, 478.69, 538.96, 618.26, 686.10, 755.74, 854.77, 918.03, 4120.8857, 4426.2296]; % (eV)
%energyAr = [15.75962, 27.62967];

% result initialization
Ar = zeros(length(energyAr)+1,M+1);
newAr = zeros(length(energyAr)+1, M+1);
Ar(1,:) = ones(1,M+1) * nppc; % one row of neutral
P = zeros(length(energyAr),M+1);
coefn = zeros(length(energyAr),1);
coef1 = zeros(length(energyAr),1);
coef2 = zeros(length(energyAr),M+1);
W = zeros(length(energyAr),M+1);

Znum = 18;
Z = 1;
UH = 13.6;
wa = 4.13e16;
Ea = 510.;
e = exp(1.);

% independent from Ea; same along entire z
% for ion = 1:length(energyAr)
%     energy = energyAr(ion);
%     coefn(ion) = 3.69*10/sqrt(energy); %n^* = 3.69*Z/sqrt(energy) 
%     coef1(ion) = 1.52e+15*(4^coefn(ion))*energy/(coefn(ion)*gamma(2*coefn(ion)));
% end
for ion = 1:length(energyAr)
    energy = energyAr(ion);
    nstar(ion)  = ion * sqrt(UH/energy);
    coef1(ion) = wa*sqrt(3)*(e/pi)^1.5 * (Z^2/nstar(ion)^4.5);
end

% insert source and update field
for n = 1:1300
    
    ndt = n*dt;
    mt = sin(2*pi*f0*ndt); % modulation
    gt = exp(-(ndt-tc)^2/(2*sigma^2)); % gaussian pulse
    source = mt*gt;
    Ex(1) = source;
    Hy = Hy-dt./mu0.*diff(Ex)/dz;
    Ex(2:M) = Ex(2:M)-dt./ep0.*diff(Hy/dz);
    
    El = abs(E0*Ex);
    
    
    % only calculate when Ea>0
    update = 0;
    for grid = 1:M+1
        if El(grid)-1e-6>0
            update = grid; % records the maximum grid number reached this condition
            newAr(:,grid) = 0;

            for ion = 1:length(energyAr)
                energy = energyAr(ion);
                coef2(ion,grid) = (4*e*Ea*ion^3/El(grid)/nstar(ion)^4)^(2*nstar(ion)-1.5);
                W(ion,grid) = coef1(ion)* coef2(ion,grid)*exp(-2*Ea/(3*El(grid)) * (ion/nstar(ion))^3);%*dt;
            end
            
            %%%% ionization solution
            % --- based on explicit solution of differential equation ---
            %%% example for 3by3 system
%             % n^0
%             newAr(1,grid) = Ar(1,grid)*exp(-W(1,grid)*dt);   
%             
%             % n^+, n^2+
%             W12 = W(1,grid)-W(2,grid);
%             if abs(W12)>eps
%                 newAr(2,grid) = W(1,grid)/W12*(exp(-W(2,grid)*dt)-exp(-W(1,grid)*dt))*Ar(1,grid) + exp(-W(2,grid)*dt)*Ar(2,grid);
%                 newAr(3,grid) = W(2,grid)/W12*exp(-W(1,grid)*dt)*Ar(1,grid) - (Ar(2,grid)+W(1,grid)/W12*Ar(1,grid))*exp(-W(2,grid)*dt)+Ar(1,grid)+Ar(2,grid)+Ar(3,grid);
%             else
%                 newAr(2,grid) = Ar(2,grid);
%                 newAr(3,grid) = Ar(3,grid);
%             end

            %%% generalization for 19by19 (or any other size)
            %%% system for the current time step

            Wmatrix = ConstructMatrix(W(:,grid));
            [Wvec, Wi] = eig(Wmatrix);
            C = Wvec\Ar(:,grid);

            for ion = 1:length(energyAr)+1
                for i = 1:length(Wi)
                    newAr(ion,grid) = newAr(ion,grid) + C(i)*Wvec(ion,i)*exp(Wi(i,i)*dt);
                end
            end
        end
    end
    
    % %%% check propagation
    %    plot(Ex)
    %    legend(int2str(n))
    %    pause(1/1000)
    %    surf(coef2)
    %        surf(exp(-W*dt))
    %        zlabel('W (s^{-1})')
    %        pause(1/1000)
    
    
    Ar(:,1:update) = newAr(:,1:update);
    %%%% check ionization
    %      surf(Ar)
    %      zlabel('no. macro particles')
    %      %pause(1/1000)
 
end

figure

subplot(3,1,1);plot(z, Ex(1:M), 'linewidth',1)
subplot(3,1,2);plot(z, W(1,1:M), 'linewidth',1)
subplot(3,1,3);plot(z, W(2,1:M), 'linewidth',1)
subplot(3,1,1);ylabel('Normalized E Field')
set(gca,'fontsize',15)
subplot(3,1,2);ylabel('W_0 (s^{-1})')
set(gca,'fontsize',15)
subplot(3,1,3);ylabel('W_1 (s^{-1})')
set(gca,'fontsize',15)
xlabel('longitudinal distance (m)','FontSize',15)

figure
subplot(4,1,1);plot(z,Ex(1:M), 'linewidth',1)
subplot(4,1,2);plot(z,Ar(1,1:M), 'linewidth',1)
subplot(4,1,3);plot(z,Ar(2,1:M), 'linewidth',1)
subplot(4,1,4);plot(z,Ar(3,1:M), 'linewidth',1)
subplot(4,1,1);ylabel('Normalized E Field','FontSize',20)
set(gca,'fontsize',15)
subplot(4,1,2);ylabel('n_0')
set(gca,'fontsize',15)
subplot(4,1,3);ylabel('n_1')
set(gca,'fontsize',15)
subplot(4,1,4);ylabel('n_2')
set(gca,'fontsize',15)
xlabel('longitudinal distance (m)','FontSize',15)
%ylabel('ionization level')

figure
plot(z,Ar(1,1:M), z,Ar(2,1:M), z,Ar(3,1:M), z,Ar(4,1:M), z,Ar(5,1:M), z,Ar(6,1:M), z,Ar(7,1:M), z,Ar(8,1:M), 'linewidth',2)
hold on
plot(z,Ex(1,1:M)*15+20,'Color',[0.25,0.25,0.25],'LineStyle',':')
set(gca,'FontSize',20);
legend('n0','n1','n2','n3','n4','n5','n6','n7','E');
xlabel('longitudinal distance (m)')



function RateMatrix = ConstructMatrix(W)
% construct Matrix of ionization rate
n = length(W)+1;
A = zeros(n,n);
W(W<eps) = 0;
% fill diagonal elements
A(1:n+1:end) = [-W;0];
% fill subdiagonal elements
A(2:n+1:end) = W;
RateMatrix = A;
end





