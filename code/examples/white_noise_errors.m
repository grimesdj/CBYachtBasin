%% White-noise record and effect of smoothing
clear; close all; clc

%% ------------------------------------------------------------
% Parameters
% -------------------------------------------------------------
N       = 2^15;                % record length
sigma_u = 3;                   % target RMS [cm/s]
win     = [2 4 8 16 32 64 128 256];

rng(1);                        % reproducible random realization


%% ------------------------------------------------------------
% 1. Construct a white spectrum with random phase
% -------------------------------------------------------------
%
% For white noise, all nonzero Fourier components have equal
% expected power. We assign equal amplitude and random phase.
%
% MATLAB FFT ordering:
%
% X(1)       = zero frequency
% X(2:N/2)   = positive frequencies
% X(N/2+1)   = Nyquist frequency
% X(N/2+2:N) = negative frequencies
%
% Hermitian symmetry is imposed so that ifft(X) is real.

X = zeros(N,1);

% Positive frequencies excluding DC and Nyquist
kpos = 2:N/2;

% Random phases
phi = 2*pi*rand(length(kpos),1);

% Unit-amplitude white spectrum
X(kpos) = exp(1i*phi);

% Hermitian-symmetric negative frequencies
X(N-kpos+2) = conj(X(kpos));

% Zero DC -> zero mean
X(1) = 0;

% Set Nyquist component to zero for convenience
X(N/2+1) = 0;


%% ------------------------------------------------------------
% 2. Inverse FFT to obtain time-domain white noise
% -------------------------------------------------------------

u = real(ifft(X));

% Remove any numerical residual mean
u = u - mean(u);

% Normalize to exactly sigma_u = 3 cm/s
u = sigma_u * u/rms(u);


%% Check statistics
fprintf('Original record\n')
fprintf('--------------------------\n')
fprintf('Mean = %8.5f cm/s\n',mean(u))
fprintf('RMS  = %8.5f cm/s\n',rms(u))
fprintf('Std  = %8.5f cm/s\n',std(u,1))
fprintf('Var  = %8.5f (cm/s)^2\n\n',var(u,1))


%% ------------------------------------------------------------
% 3. Estimate power spectrum
% -------------------------------------------------------------
%
% With sample spacing dt = 1, frequency is expressed in
% cycles/sample. Change dt below if a physical sample rate is known.

dt = 1;
fs = 1/dt;

U = fft(u);

% Two-sided PSD
Suu2 = (abs(U).^2)/(fs*N);

% One-sided PSD
Suu = Suu2(1:N/2+1);
Suu(2:end-1) = 2*Suu(2:end-1);

f = (0:N/2)'*fs/N;


%% ------------------------------------------------------------
% 4. Smooth the record with different boxcar lengths
% -------------------------------------------------------------

AVG = zeros(size(win));
RMS = zeros(size(win));

for ii = 1:length(win)

    nw = win(ii);

    % Centered moving average
    us = movmean(u,nw);

    % Statistics of smoothed record
    AVG(ii) = mean(us);
    RMS(ii) = rms(us);

end


%% ------------------------------------------------------------
% Theoretical RMS for averaging independent white-noise samples
% -------------------------------------------------------------
RMS_theory = sigma_u ./ sqrt(win);


%% ------------------------------------------------------------
% Display results
% -------------------------------------------------------------
T = table(win(:),AVG(:),RMS(:),RMS_theory(:), ...
    'VariableNames',{'Window','AVG_cm_s','RMS_cm_s','RMS_theory_cm_s'});
disp(T)


%% ------------------------------------------------------------
% 5. Plot example white-noise record
% -------------------------------------------------------------
figure
nplot = 1:500;
plot(nplot,u(nplot),'k-')

xlabel('Sample')
ylabel('u (cm s$^{-1}$)')
title('White-noise velocity record')
grid on


%% ------------------------------------------------------------
% 6. Plot power spectrum
% -------------------------------------------------------------
figure
loglog(f(2:end),Suu(2:end),'k-')

xlabel('Frequency (cycles/sample)')
ylabel('PSD [(cm s$^{-1}$)$^2$/(cycles/sample)]')
title('White-noise power spectrum')
grid on


%% ------------------------------------------------------------
% 7. Plot AVG and RMS versus smoothing-window length
% -------------------------------------------------------------

figure
yyaxis left

loglog(win,abs(AVG),'o-','LineWidth',1.5,'MarkerSize',7)
ylabel('AVG (cm s$^{-1}$)')
ax = gca;
ax.YAxis(1).Scale='log';

yyaxis right

loglog(win,RMS,'s-','LineWidth',1.5,'MarkerSize',7)
hold on

loglog(win,RMS_theory,'k--','LineWidth',1.5)

ylabel('RMS (cm s$^{-1}$)')

xlabel('Smoothing window, N_w (samples)')

legend('AVG','RMS','3/\surdN_w', ...
       'Location','best')

grid on
box on

title('Effect of smoothing white noise')