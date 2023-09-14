% This script implements the S-DSM controller with dynamic alpha proposed in the paper 
% "A Stochastic Controller for Primary Frequency Regulation using ON/OFF Demand Side Resources"
% Guanyu Tian 03/07/2023

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
WH_T_record = zeros(N,N_step);% random initial temeprature
WH_T_record(:,1) = rand(N,1)*10+50;
% WH_T_record(1:50,1) = ones(50,1) * (50 + 0.02);
% WH_T_record(N-49:N,1) = ones(50,1) * (60 - 0.02);
Pe_fix = Pm-sum(WH_status*P_on)*ones(1,N_step);
Pe_fix_init = Pe_fix(1);
k=10;

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
    alpha = abs(omega-1)*k;
     
    if mod(t_sequence(t),0.1)==0 % response frequency of WH is 0.1 s
        if omega<oemga_DR(1) % WH response to under frequency
            for i = 1:N
                if WH_status(i)==1
                    n = rand(1);
                    if n<alpha
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
    
    WH_T_record(:,t) = WH_T_record(:,t-1) + 0.01/120*(2*WH_status-1);
    
    % temperature trigger
    index_ON = find(WH_T_record(:,t)<50);
    index_OFF = find(WH_T_record(:,t)>60);
    WH_status(index_ON) =  WH_status(index_ON)*0+1;
    WH_status(index_OFF) =  WH_status(index_OFF)*0;
    
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
legend('Benchmark','S-DSM-static \alpha','S-DSM-dynamic \alpha')
grid on
title(['#GIWH=' num2str(N)])
subplot(2,1,2)
hold on
plot(t_sequence,sum(WH_status_record)/N,'LineWidth',3)
xlabel('Time(s)')
ylabel('Ratio of WH is ON')
legend('Benchmark','S-DSM-static \alpha','S-DSM-dynamic \alpha')
grid on
% ylim([0 0.6])
end

% count switching times for N=2800
sw_count = zeros(N,1);
for m_index = 1:N
    for t = 2:N_step
        if WH_status_record(m_index,t) ~= WH_status_record(m_index,t-1)
            sw_count(m_index) = sw_count(m_index) + 1;
        end
    end
end
prob = [length(find(sw_count==0))/N;
        length(find(sw_count==1))/N;
        length(find(sw_count==2))/N;
        length(find(sw_count==3))/N;
        length(find(sw_count==4))/N;
        length(find(sw_count==5))/N]*100;
list = 0:3;
figure
hold on
for i=list
    line([i i],[0 prob(i+1)],'LineWidth',2,'Color','black')
end
xlim([-0.5 max(list)+0.5])
ylim([0 max(prob+10)])
xlabel('Number of Switching')
ylabel('Ratio (%)')
scatter(list,prob(list+1),'MarkerEdgeColor','red',...
              'MarkerFaceColor','red',...
              'LineWidth',2)
grid on
title('Distribution of number-of-switching')

























