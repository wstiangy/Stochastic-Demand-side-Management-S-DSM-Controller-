% This script implements the deterministic controller with dead bands for the paper 
% "A Stochastic Controller for Primary Frequency Regulation using ON/OFF Demand Side Resources"
% Guanyu Tian 02/28/2023

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
AR = 2.5*[ones(1,2/0.01+1)*0.2 ones(1,2/0.01)*0.1 ones(1,6/0.01)*0.03;
      -ones(1,2/0.01+1)*0.2 -ones(1,2/0.01)*0.1 -ones(1,6/0.01)*0.03];% define dead bands
% figure
% hold on
% plot(t_sequence, 60+AR(1,:))
% plot(t_sequence, 60+AR(2,:))
N_step = length(t_sequence);
WH_status = [zeros(N/2,1);ones(N/2,1)]; % half on, half off
WH_status_record = repmat(WH_status,[1,N_step]);
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
    oemga_DR = omega0-AR(:,t)/60;
    
    if mod(t_sequence(t),0.1)==0 % response frequency of WH is 0.1 s
        delay_index = min(t + randi([0 1000],1,N),N_step); % randomize response time
%         delay_index = ones(1,N)*t; % instant response
        if omega<oemga_DR(1) % WH response to under frequency
%             WH_status = zeros(N,1); % turn off all WH immediately
            % add random delay into response
            for i = 1:N
                WH_status_record(i,delay_index(i):N_step) = zeros(1,N_step-delay_index(i)+1);
            end
        end
        if omega>oemga_DR(2) % WH response to under frequency
%             WH_status = ones(N,1); % turn on all WH immediately
            % add random delay into response
            for i = 1:N
                WH_status_record(i,delay_index(i):N_step) = ones(1,N_step-delay_index(i)+1);
            end
        end
    end
    Pe = Pe_fix(t) + sum(WH_status_record(:,t)*P_on); 
%     WH_status_record(:,t) = WH_status;
        
    ddelta = omega-omega0;
    domega = (Pm-Pe-D*(omega-omega0))/(2*H);
    
    x = [x [delta+ddelta;omega+domega]];
end

figure(N)
subplot(2,1,1)
hold on
plot(t_sequence(1:end-1),x(2,1:end-1)*60,'LineWidth',3)
xlabel('Time(s)')
ylabel('Frequency(Hz)')
% legend('N=2800','N=2200','N=1800')
grid on
subplot(2,1,2)
hold on
plot(t_sequence(1:end-1),sum(WH_status_record(:,1:end-1))/N,'LineWidth',3)
xlabel('Time(s)')
ylabel('Ratio of WH is ON')
% legend('N=2800','N=2200','N=1800')
grid on
% ylim([0 0.6])

end

























