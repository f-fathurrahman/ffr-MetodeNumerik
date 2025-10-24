clc; close all; clear all;
% Construct a Gaussian Waveform in time, frequency spectrum of which 
% has its magnitude at 1 GHz as 10% of the maximum.
maximum_frequency = 1e9;
tau = sqrt(2.3)/(pi*maximum_frequency);
t_0 = 4.5 * tau;
time_array = [1:1000]*1e-11;
g = exp(-((time_array - t_0)/tau).^2);
figure(1); 
plot(time_array*1e9, g,'b-','linewidth',1.5);
title('g(t)=e^{-((t-t_0)/\tau)^2}','fontsize',14);
xlabel('time (ns)','fontsize',12);
ylabel('magnitude','fontsize',12);
set(gca,'fontsize',12);
grid on;

% Perform time to frequency domain transform 
frequency_array = [0:1000]*2e6;
dt = time_array(2)-time_array(1);
G = time_to_frequency_domain(g, dt, frequency_array, 0);

figure(2); 
subplot(2,1,1);
plot(frequency_array*1e-9, abs(G),'b-','linewidth',1.5);
title('G(\omega) = F(g(t))','fontsize',12);
xlabel('frequency (GHz)','fontsize',12);
ylabel('magnitude','fontsize',12);
set(gca,'fontsize',12);
grid on;
subplot(2,1,2);
plot(frequency_array*1e-9, angle(G)*180/pi,'r-','linewidth',1.5);
xlabel('frequency (GHz)','fontsize',12);
ylabel('phase (degrees)','fontsize',12);
set(gca,'fontsize',12);
grid on;
drawnow;

% Perform frequency to time domain transform 
df = frequency_array(2)-frequency_array(1);
g2 = frequency_to_time_domain(G, df, time_array);

figure(3); 
plot(time_array*1e9, abs(g2),'b-','linewidth',1.5);
title('g(t)= F^{-1}(G(\omega))','fontsize',14);
xlabel('time (ns)','fontsize',12);
ylabel('magnitude','fontsize',12);
set(gca,'fontsize',12);
grid on;

