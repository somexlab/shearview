clear; close all; clc;


path = input('Enter path:  '); % enter path of the data file

A = readmatrix(path);
A(end,:) = [];
% write results to txt file
file = 'your_file.txt'; % file where the resluts will be saved


%--------Parameters-------------%

volt=0.689;%voltage to distance constant (V/mm)
Fconstant = 0.911;% Force constant of voice coil (N/A)
fs = 500; % Sampling frequency (Hz) 
targetSampleRate = 500;
gap=0.65*1e-3; %m
area=560*1e-6; %m^2, calculated from the cross section of the sample with a picture taken with a smartphone
% passband = (omega/2*pi) / 20;

%--------Parameters-------------%


time=(A(:,2));
% make sure the voltage sine is correct
X1=(A(:,3)); %raw displacement data in Volts unit
current=(A(:,7));
% T = 2; % Total duration of simulated signal (s)

Dx=X1; 
Dx=Dx-Dx(1,1);

Dx=(Dx/volt)*1e-3;%displacement in m 
% strain=Dx/gap;% (-)
figure
plot(time,Dx,LineWidth=2)

%select the regime of data to continue with the analysis.

[x,y]=ginput(2);
%get the indices of the selected regime
v1=x(1,1); v2=x(2,1);
%find the part of the array in time (never noisy) that is within the regime
[r1,c1]=find(time>=v1); [r2,c2]=find(time<=v2);

%get the data to continue 
time=time(r1(1,1):r2(end,1),1);
time=time-time(1,1);
Dx=Dx(r1(1,1):r2(end,1),1);
current=current(r1(1,1):r2(end,1),1);

%substract DC components if present
Dx = Dx - mean(Dx);
current = current - mean(current);

% --- Plot Raw Data ---
figure
yyaxis left
plot(time,Dx,LineWidth=1.2)
ylabel('Displacement (m)',FontSize=15)
yyaxis right
plot(time,current,LineWidth=1.2)
ylabel('Current (A)',FontSize=15)
xlabel('Time (s)',FontSize=15)
title('raw data',FontSize=15)

% Signals do not have a uniform sampling rate due to the control loop submillisecond 
% timing errors. Resample the two signals to a common timebase


[rDx,trDx] = resample(Dx,time,targetSampleRate,'spline');
[rcurrent, trcurrent] = resample(current,time,targetSampleRate,'spline');
rtime=trDx;

% Apply windowing to the signals
window_type = 'hann';
window_x = feval(window_type, length(rDx));
window_I = feval(window_type, length(rcurrent));

x_windowed = rDx .* window_x;
I_windowed = rcurrent .* window_I;

% Calculate window gain for amplitude correction
% Coherent gain (sum of window values / N)

window_gain_x = sum(window_x) / length(x_windowed);
window_gain_I = sum(window_I) / length(I_windowed);

figure;
subplot(2,1,1);
plot(trDx, x_windowed, 'b', 'DisplayName', 'Windowed strain(t)'); hold on;
xlabel('Time (s)');
ylabel('Amplitude');
title('Strain Signal Windowing');
legend show;
grid on;
hold off;

subplot(2,1,2);
plot(trcurrent, I_windowed, 'r', 'DisplayName', 'Windowed stress(t)');
xlabel('Time (s)');
ylabel('Amplitude');
title('Stress Signal Windowing');
legend show;
grid on;
hold off;

% FFT Analysis for Amplitude and Phase

N = length(x_windowed); % Number of samples
L = trDx(end,1); % Length of signal (s) - assuming T is total duration
f = fs*(0:(N/2))/N; % Frequency vector for single-sided spectrum

% Perform FFT
X_fft = fft(x_windowed);
I_fft = fft(I_windowed);

% Compute the single-sided spectrum

P1_x = X_fft(1:N/2+1);
P1_x(2:end-1) = 2*P1_x(2:end-1); % Multiply by 2 for non-DC, non-Nyquist components
% Apply window gain correction for amplitude
P1_x = P1_x / window_gain_x; % Correct for amplitude reduction due to window
% P1_power = abs(X_fft).^2 / window_gain_x; % single-sided power spectrum
% figure
% P2 = abs(Y).^2;
% P1 = P2(1:N/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% plot(f,P1_power(0:(N/2))/N)

P1_I = I_fft(1:N/2+1);
P1_I(2:end-1) = 2*P1_I(2:end-1); % Multiply by 2 for non-DC, non-Nyquist components
P1_I = P1_I / window_gain_I; % Correct for amplitude reduction due to window

% --- Find Dominant Frequency ---
[~, idx_dominant_x] = max(abs(P1_x(2:end))); % Exclude DC component
idx_dominant_x = idx_dominant_x + 1; % Adjust index since we started from 2
f_dominant_x = f(idx_dominant_x);

[~, idx_dominant_I] = max(abs(P1_I(2:end))); % Exclude DC component
idx_dominant_I = idx_dominant_I + 1; % Adjust index since we started from 2
f_dominant_I = f(idx_dominant_I);

fprintf('\n--- FFT Analysis Results ---\n');
fprintf('Dominant Frequency (Displacement): %.2f Hz\n', f_dominant_x);
fprintf('Dominant Frequency (Current): %.2f Hz\n', f_dominant_I);

% --- Calculate Amplitude ---
A_calculated = abs(P1_x(idx_dominant_x)) / N; % Magnitude from FFT is sum of sines, /N normalizes to average
fprintf('Calculated Amplitude Strain (FFT): %.5f\n', A_calculated);

I_calculated = abs(P1_I(idx_dominant_I)) / N; % Magnitude from FFT is sum of sines, /N normalizes to average
fprintf('Calculated Amplitude Current(FFT): %.5f\n', I_calculated);
% --- Calculate Phase ---
% Phase is given by the angle of the complex FFT value at the dominant frequency.
% MATLAB's angle() function returns radians.
% The phase from FFT corresponds to a cosine term (A*cos(wt + phi_cos)).
% If your model is A*sin(wt + phi_sin), then phi_sin = phi_cos - pi/2.

phi_calculated_cos_rad = angle(X_fft(idx_dominant_x));
phi_calculated_sin_rad = phi_calculated_cos_rad - pi/2;

fprintf('Calculated Phase (relative to cos, FFT): %.2f rad (%.2f degrees)\n', ...
    phi_calculated_cos_rad, rad2deg(phi_calculated_cos_rad));
fprintf('Calculated Phase (relative to sin, FFT): %.2f rad (%.2f degrees) (True: %.2f degrees)\n', ...
    phi_calculated_sin_rad, rad2deg(phi_calculated_sin_rad));

% --- Calculate Phase Difference between Displacement and Current ---
phi_X_rad = angle(X_fft(idx_dominant_x)); % Phase of displacement at dominant freq
phi_I_rad = angle(I_fft(idx_dominant_I)); % Phase of current at dominant freq

phase_diff_rad = phi_I_rad - phi_X_rad ;
% Ensure phase difference is wrapped to [-pi, pi]
phase_diff_rad = atan2(sin(phase_diff_rad), cos(phase_diff_rad)); % Wraps to [-pi, pi]

fprintf('Phase Difference (Stress rel. to Strain): %.2f rad (%.2f degrees)\n', ...
    phase_diff_rad, rad2deg(phase_diff_rad));
% fprintf(' (Positive means displacement leads current, negative means displacement lags current)\n');
% fprintf(' (True phase diff: %.2f degrees, calculated as (pi/4 - pi/8) = pi/8 = 22.5 deg)\n', rad2deg(phi_true_rad - I_phase_rad));


% Find the index of the third harmonic for strain
f_third_strain = 3 * f_dominant_x;
[~, idx_third_x] = min(abs(f - f_third_strain)); % Find frequency closest to 3 * dominant freq
% Ensure it's not the dominant frequency itself if there's an exact match issue
if idx_third_x == idx_dominant_x && length(f) > idx_third_x
    idx_third_x = idx_third_x + 1; % Move to next frequency bin if it coincides
end

% Find the index of the third harmonic for stress
f_third_stress = 3 * f_dominant_I;
[~, idx_third_I] = min(abs(f - f_third_stress)); % Find frequency closest to 3 * dominant freq
% Ensure it's not the dominant frequency itself
if idx_third_I == idx_dominant_I && length(f) > idx_third_I
    idx_third_I = idx_third_I + 1; % Move to next frequency bin
end

% Get amplitudes of the third harmonics
A_third_strain = abs(P1_x(idx_third_x)) / N;
A_third_stress = abs(P1_I(idx_third_I)) / N;

fprintf('Third Harmonic Frequency (Strain): %.2f Hz (Index: %d)\n', f(idx_third_x), idx_third_x);
fprintf('Third Harmonic Amplitude (Strain): %.5f\n', A_third_strain);
fprintf('Third Harmonic Frequency (Stress): %.2f Hz (Index: %d)\n', f(idx_third_I), idx_third_I);
fprintf('Third Harmonic Amplitude (Stress): %.5f\n', A_third_stress);


% Calculate the ratios
    ratio_strain = A_third_strain / A_calculated;
    fprintf('Ratio (Strain 3rd/1st Harmonic): %.4f\n', ratio_strain);

    ratio_stress = A_third_stress / I_calculated;
    fprintf('Ratio (Stress 3rd/1st Harmonic): %.4f\n', ratio_stress);

% --- Plot Magnitude Spectrum ---
figure;
subplot(2,1,1);
plot(f, abs(P1_x)/N); % Plot normalized amplitude spectrum
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Strain Magnitude Spectrum (FFT)');
% xlim([0 f_true*2.5]); % Limit x-axis to relevant frequencies
grid on;

subplot(2,1,2);
plot(f, abs(P1_I)/N); % Plot normalized amplitude spectrum
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Stress Magnitude Spectrum (FFT)');
% xlim([0 f_true*2.5]); % Limit x-axis to relevant frequencies
grid on;

% --- Plot Phase Spectrum ---
figure;
subplot(2,1,1);
plot(f, rad2deg(angle(P1_x))); % Plot phase in degrees
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Strain Phase Spectrum (FFT)');
% xlim([0 f_true*2.5]);
ylim([-180 180]); % Phase is typically plotted from -180 to 180
grid on;

subplot(2,1,2);
plot(f, rad2deg(angle(P1_I))); % Plot phase in degrees
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Stress Phase Spectrum (FFT)');
% xlim([0 f_true*2.5]);
ylim([-180 180]);
grid on;

% Calculate the power spectrum of stress
P2 = abs(I_fft).^2;
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
P1(P1==0) = eps;
P1_dB = 10*log10(P1);
%%
figure;
plot(f, P1_dB);
title('Power Spectrum of Stress');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
grid on;
xlim([-10 100])
%%

%------ Calculation of moduli---------------------
mass=0.246; % kg (measured the moving parts with a scale)
omega=2*pi*f_dominant_x;
I=mass*(gap/area);

stress_amp=(I_calculated*Fconstant)/(area); % Pa
strain_amp=A_calculated/gap;%the amplitude of the fitting of strain data

G1_corr=I*omega^2; 
G1=(stress_amp/strain_amp)*cos(phase_diff_rad);
Gprime=G1+G1_corr;
% F_fr0=0.14*1e-3; %(N) sliding friction
% G2_corr=F_fr0/(strain_amp*area);
G2=(stress_amp/strain_amp)*sin(phase_diff_rad);
G2prime=G2;
%%
% write data in a file
T=['Angular_Frequency(rad/s) strain_amp(-) stress_amp(Pa) current(A) phase_shift_angle(rad) Gprime(Pa) G2prime(Pa)'];
formatSpec = '%.5e %.5e %.5e %.5e %.3e %.3e %.3e %.3e %.3e\n';
fileID = fopen(file,"a+");
fprintf(fileID,formatSpec,T);
fclose(fileID);


results = [omega strain_amp stress_amp I_calculated phase_diff_rad Gprime G2prime ratio_strain ratio_stress];
 

% Open once at the beginning in append mode
fid = fopen(file, 'a');

% Optional: write header only if file is empty
if ftell(fid) == 0
    fprintf(fid, 'Angular_Frequency(rad/s)\tstrainAmplitude(-)\tstressAmplitude(Pa)\tcurrentAmplitude(A)\tphase_shift_angle(rad)\tGprime(Pa)\tG2prime(Pa)\tratioStrain\tratioStress\n');
end

    fprintf(fid, '%.5f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\t%.7f\n', results);

% Close the file when done
fclose(fid);



