%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECEPTION USING A LINEAR ANTENNA ARRAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ARRAY PARAMETERS
c = 3e8; % propagation velocity  [m/s]
f0 = 1e9; % carrier frequency [Hz]
lambda = c/f0; % wavelength [m]

La = 3;         % Array aperture [m]
N_ant = 30;     % Number of elements
delta_z = La/N_ant; % spacing between neighbouring elements
zn = (-(N_ant-1)/2:(N_ant-1)/2)*delta_z; % vertical position of each element

fprintf("Minimum number of elements required for %.2fm array (180deg FOV, avoid replicas): %d\n  Selected: %d antennas\n", La, ceil(2*La/lambda), N_ant);
fprintf("Minimum angle resolution with this configuration: %.2fdeg\n", rad2deg(lambda/La));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIRECTION SCANNING
% Directions tested by setting the array coefficients
theta = linspace(-90,90,1001);
% Corresponding spatial frequencies
fz = sind(theta)/lambda;
% Coefficients
an = exp(-1i*2*pi*zn(:)*fz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RECEPTION OF A SHORT PULSE IN THE TIME DOMAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PULSE PARAMETERS -note: see function rec_sign to change the type of signal
B = 1e3; % Bandwidth [Hz]
T_sig = 1/B;
T_obs = 30*T_sig; % Total observation time [s]
dt = 1/4*T_sig; % sampling time
t = (-T_obs/2:dt:T_obs/2); % time axis

% Incident direction
theta_1 = 0;% [deg]
theta_2 = 5.73; % [deg]
%theta_3 = 60;   % if you want to add a third satellite

% Delays with respect to the antenna in z = 0
tau = 1/c * sind(theta_1)*zn(:);
tau2 = 1/c*sind(theta_2)*zn(:); % delay of second satellite
%tau3 = 1/c*sind(theta_3)*zn(:); % delay of third satellite
deltaT = 3e-3; % [s] time delay between one signal to another (by default the first signal is always in 0)

% plotting both signals
plotsignal(t,rec_sign(t, T_sig, 0));
plotsignal(t,rec_sign(t, T_sig, 0));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Signals received at the N antennas, for each time
sn = zeros(N_ant,length(t));

for ant_id = 1:N_ant % for each antenna, this computes 
    sn(ant_id,:) = rec_sign(t, T_sig, tau(ant_id)) * exp(-1i*2*pi*f0*tau(ant_id));
    sn(ant_id,:) = sn(ant_id,:) + rec_sign(t, T_sig, tau2(ant_id)) * exp(-1i*2*pi*f0*tau2(ant_id));
    %sn(ant_id,:) = sn(ant_id,:) + rec_sign(t, T_sig, tau3(ant_id)+deltaT)* exp(-1i*2*pi*f0*tau3(ant_id)); % third sat
end
% !!!!!! NOTE THAT: f0*tau(n) = 1/lambda*sind(theta_i)*zn(:);
% => SAME PHASE AS IN THE FREQUENCY DOMAIN APPROACH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(2,1,1), imagesc(zn,t,abs(sn')), axis xy,colorbar
hold on
plot(zn,min(t)*ones(size(zn)),'rv','LineWidth',3)
ylabel('time [s]'), xlabel('antenna position [m]')
title('Received signals (abs)')
subplot(2,1,2), imagesc(zn,t,angle(sn')), axis xy,colorbar
hold on
plot(zn,min(t)*ones(size(zn)),'rv','LineWidth',3)
ylabel('time [s]'), xlabel('antenna position [m]')
title('Received signals (phase)')

% Combined signal (for each time instant) - compute Fourier transform
S = an'*sn;
% remember:
%    an : let's say it is the "spatial delay" - it represents the sampling
%         distance/position of each element
%    sn : vector of signals received at each receiver
%
%     S : Discrete Fourier transform of the signal, remember that:
%           F(sn) = sum( sn * exp(-1i * 2*pi * fz * dz) )
%          such to obtain N_ant elements for S
figure
imagesc(theta,t,abs(S')), axis xy,colorbar
hold on
ylabel('time [s]'), xlabel('\theta [deg]'), 
title('Array response')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% array response at t = 0s; - plot the fourier transform received at t=0
figure
plot(theta,abs(S(:,t == 0)'), 'LineWidth', 1.5)
hold on, grid on, grid minor
ylabel('signal amplitude'), xlabel('\theta [deg]'), 
title('Array response for t = 0s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% introducing a function to change the received signal
function signal = rec_sign(time, T_sig, del)
    % signal = rec_sign(time, T_sig, del)
    %  requires:
    %    - time  : time array to evaluate signal on
    %    - T_sig : T_sig = 1/Bandwidth
    %    - del   : delay to apply
    %

    s_sel = 3; % select the signal
    
    switch s_sel
        case 1
            signal = rectpuls((time-del)/T_sig) + noise(time,.1);   
        case 2
            signal = rectpuls((time-del)/T_sig) + noise(time,.1); % using a rectangular pulse
        case 3
            signal = sinc((time-del)/T_sig);% + noise(time,.2);       % Tebaldini's default: the received signal is a cardinal sine, delayed of tau.
        otherwise
            error("Selected Signal Unavailable");
    end
end