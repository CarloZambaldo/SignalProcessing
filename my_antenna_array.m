%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECEPTION USING A LINEAR ANTENNA ARRAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ARRAY PARAMETERS
c = 3e8; % propagation velocity  [m/s]
f0 = 1.2e9; % carrier frequency [Hz]
lambda = c/f0; % wavelength [m]

La = 3; % Array aperture [m]
N_ant = 30; % Number of elements
fprintf("Minimum number of elements required for %.2fm array (180deg FOV, avoid replicas): %d\n  Selected: %d\n", La, ceil(2*La/lambda), N_ant);
fprintf("Minimum angle resolution with this configuration: %.2fdeg\n", rad2deg(lambda/La));
delta_z = La/N_ant; % spacing between neighbouring elements
zn = (-(N_ant-1)/2:(N_ant-1)/2)*delta_z; % vertical position of each element

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIRECTION SCANNING
% Directions tested by setting the array coefficients
theta = linspace(-90,90,501);
% Corresponding spatial frequencies
fz = sind(theta)/lambda;
% Coefficients
an = exp(-1i*2*pi*zn(:)*fz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECEPTION OF A SHORT PULSE IN THE TIME DOMAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PULSE PARAMETERS
B = 1e6; % Bandwidth [Hz]
T_sig = 1/B;
T_obs = 20*T_sig; % Total observation time [s]
dt = 1/4*T_sig; % sampling time
t = (-T_obs/2:dt:T_obs/2); % time axis

% Incident direction
theta_1 = 30+3;    % [deg]
theta_2 = 30+-6; % [deg]

% Delays with respect to the antenna in z = 0
tau = 1/c * sind(theta_1)*zn(:);
tau2 = 1/c*sind(theta_2)*zn(:); % delay of second satellite, imposing the second sat at 30°
deltaT = 0;1e-3; % [s] time delay between one signal to another (by default the first signal is always in 0)

% plotting both signals
plotsignal(t,rec_sign(t, T_sig, 0));
plotsignal(t,rec_sign(t, T_sig, deltaT));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signals received at the N antennas, for each time
sn = zeros(N_ant,length(t));

for ant_id = 1:N_ant % for each antenna, this computes 
    sn(ant_id,:) = rec_sign(t, T_sig, tau(ant_id)) * exp(-1i*2*pi*f0*tau(ant_id));
    sn(ant_id,:) = sn(ant_id,:) + rec_sign(t, T_sig, tau2(ant_id)+deltaT) * exp(-1i*2*pi*f0*tau2(ant_id));
end
% !!!!!! NOTE THAT: f0*tau(n) = 1/lambda*sind(theta_i)*zn(:);
% => SAME PHASE AS IN THE FREQUENCY DOMAIN APPROACH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(2,1,1), imagesc(zn,t,abs(sn')), axis xy,colorbar
hold on
plot(zn,min(t)*ones(size(zn)),'ro','LineWidth',5)
ylabel('time [s]'), xlabel('antenna position [m]')
title('Received signals (abs)')
subplot(2,1,2), imagesc(zn,t,angle(sn')), axis xy,colorbar
hold on
plot(zn,min(t)*ones(size(zn)),'ro','LineWidth',5)
ylabel('time [s]'), xlabel('antenna position [m]')
title('Received signals (phase)')

% Combined signal (for each time instant) - compute Fourier transform
S = an'*sn;
% remember:
%    an : let's say it is the "spatial delay" - it represents the sampling
%         distance/position of each element
%    sn : vector of signals received at each receiver
%
%     S : Fourier transform of the signal, remember that:
%           F(sn) = sum( sn * exp(-1i * 2*pi * fz * dz) )
%          such to obtain N_ant elements for S
figure
imagesc(theta,t,abs(S')), axis xy,colorbar
hold on
ylabel('time [s]'), xlabel('\theta [deg]'), 
title('Array response')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% array responce at t = 0s;

figure
plot(theta,abs(S(:,t == 0)'), 'LineWidth', 1.5)
hold on, grid on, grid minor
ylabel('signal amplitude'), xlabel('\theta [deg]'), 
title('Array response for t = 0s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

function signal = rec_sign(time, T_sig, del)
    % introducing a function to change the received signal
        %signal = rectpuls((time-del)/T_sig) + noise(time,.1);   
        % signal = rectpuls((time-del)/T_sig) + noise(time,.1); % using a rectangular pulse
        signal = sinc((time-del)/T_sig); + noise(time,.2);       % Tebaldini's default: the received signal is a cardinal sine, delayed of tau.
end