% This script implements the S-DSM controller with fixed alpha proposed in the paper 
% "A Stochastic Controller for Primary Frequency Regulation using ON/OFF Demand Side Resources"
% Guanyu Tian 09/09/2022

clear all
% close all
clc

for N = [1800, 2100, 2800] % number of water heater participating in DR, 3 scenarios
% initialization
P_on = 1e-4; % power of water heater is 1 kW when turned on
omega0 = 1; % nominal frequency is 60 Hz
oemga_DR = [60,60]/60; % the thresholds for DR triggering, 2 Hz deadband
H = 20;
D = 3;
Pm = 10;
T = 10;
t_sequence = 0:0.01:T;
N_step = length(t_sequence);
WH_status = [zeros(N/2,1);ones(N/2,1)]; % half on, half off
WH_status_record = zeros(N,N_step);
WH_status_record(:,1) = WH_status;
Pe_fix = Pm-sum(WH_status*P_on)*ones(1,N_step);
Pe_fix_init = Pe_fix(1);

% time domain simulation
rng(666)
x = [0;omega0];
flag = zeros(N,1);
for t=2:N_step
    delta = x(1,end);
    omega = x(2,end);
    if t_sequence(t)>=1
        Pe_fix(t)=Pe_fix_init-0.1 + normrnd(0,0.005);
    end
    alpha = 0.1;
     
    if mod(t_sequence(t),0.1)==0 % response frequency of WH is 0.1 s
        if omega<oemga_DR(1) % WH response to under frequency
            for i = 1:N
                if WH_status(i)==1
                    n = rand(1);
                    if n<=alpha
                        WH_status(i) = 0;
                    end
                end
            end
        end
        if omega>oemga_DR(2) % WH response to over frequency
            for i = 1:N 
                if WH_status(i)==0
                    n = rand(1);
                    if n<=alpha
                        WH_status(i) = 1;
                    end
                end
            end
        end
    end
    Pe = Pe_fix(t) + sum(WH_status*P_on); 
    WH_status_record(:,t) = WH_status;
        
    ddelta = omega-omega0;
    domega = (Pm-Pe-D*(omega-omega0))/(2*H);
    
    x = [x [delta+ddelta;omega+domega]];
end

figure(N)
subplot(2,1,1)
hold on
plot(t_sequence,x(2,:)*60,'LineWidth',3)
xlabel('Time(s)')
ylabel('Frequency(Hz)')
% legend('N=1900','N=2500','N=3100')
grid on
subplot(2,1,2)
hold on
plot(t_sequence,sum(WH_status_record)/N,'LineWidth',3)
xlabel('Time(s)')
ylabel('Ratio of WH is ON')
% legend('N=1900','N=2500','N=3100')
grid on
% ylim([0 0.6])

end


























