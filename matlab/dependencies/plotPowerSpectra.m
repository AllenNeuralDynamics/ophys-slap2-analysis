function plotPowerSpectra(VV, Fs)
% function to plot normalized power spectra of several time series
% courtesy of ChatGPT

% Assume VV is your TxN matrix where T is the number of time points
% and N is the number of time series.
[T, N] = size(VV); % Get the dimensions of the matrix
% Frequency vector (assuming sampling rate is 1, adjust as necessary)
if nargin < 2
    Fs = T/120; % Sampling frequency
end
f = (0:T/2-1)*(Fs/T); % Frequency vector
% Pre-allocate space for power spectra
powerSpectra = zeros(length(f), N);
% Prepare a figure outside the loop for all plots
figure;
hold on; % This allows all plots to be drawn on the same figure
% Loop through each time series
colors = distinguishable_colors(N);
for i = 1:N
    % Compute the Fourier Transform of the time series
    Y = fft(VV(:,i));
    % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
    P2 = abs(Y/T);
    P1 = P2(1:T/2);
    P1(2:end-1) = 2*P1(2:end-1);
    % Normalize the power spectrum
    P1_normalized = P1 / max(P1);
    % Store the normalized power spectrum in the matrix
    powerSpectra(:,i) = P1_normalized;
    % powerSpectra(:,i) = P1;
    % Plot the normalized power spectrum
    plot(f, P1_normalized, 'DisplayName', ['Trace ' num2str(i)], 'Color',[colors(i,:), 0.5]);
end
hold off; % No more plots to add
title('Power Spectra of Time Series');
xlabel('Frequency (Hz)');
ylabel('Power');
legend show; % Display a legend to identify each time
end