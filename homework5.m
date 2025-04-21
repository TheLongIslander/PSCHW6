% Homework 5 Spring 2025 - Full Simulation
% RUID: 208001821

clc; clear; close all;

% ----- Parameters -----
RUID = 208001821;
rng(RUID);

T = 2;                  % Bit duration (s)
A = 1;                  % Pulse height (Volt)
Ts = 0.02;              % Sampling interval (s)
R = 1 / T;              % Bit rate (Hz)
N = 100000;             % Number of bits
bit_samples = round(T / Ts);
L = N * bit_samples;    % Signal length
t_sinc = -5*T:Ts:5*T;   % Truncated sinc window

% ----- Generate random bits -----
bb = randi([0, 1], 1, N);
disp('First 10 bits:');
disp(bb(1:10));

% ----- Define pulses -----
p = A * ones(1, bit_samples);            % Square pulse
ps = A * sinc(t_sinc / T);               % Truncated sinc pulse
Ep = sum(p.^2) * Ts;
Eps = sum(ps.^2) * Ts;

% ----- Modulate signal -----
s = zeros(1, L + length(ps));            % Using p(t)
ss = zeros(1, L + length(ps));           % Using ps(t)
for i = 1:N
    idx = (i-1)*bit_samples + 1;
    if bb(i) == 1
        s(idx:idx+bit_samples-1) = s(idx:idx+bit_samples-1) + p;
        ss(idx:idx+length(ps)-1) = ss(idx:idx+length(ps)-1) + ps;
    else
        s(idx:idx+bit_samples-1) = s(idx:idx+bit_samples-1) - p;
        ss(idx:idx+length(ps)-1) = ss(idx:idx+length(ps)-1) - ps;
    end
end

% ----- Plot pulses -----
figure;
subplot(2,1,1); plot(p); title('Square Pulse p(t)');
subplot(2,1,2); plot(ps); title('Sinc Pulse ps(t)');

fprintf('Energy of p(t): %.4f, Energy of ps(t): %.4f\n', Ep, Eps);

% ----- Generate noise -----
SNR_dB = 0:1:7;
SNR = 10.^(SNR_dB / 10);
ni = cell(1, 8); nis = cell(1, 8);
for j = 1:8
    vi = Ep / SNR(j);
    vis = Eps / SNR(j);
    ni{j} = sqrt(vi) * randn(1, length(s));
    nis{j} = sqrt(vis) * randn(1, length(ss));
end

% ----- Add noise -----
d = cell(1, 8); ds = cell(1, 8);
for j = 1:8
    d{j} = s + ni{j};
    ds{j} = ss + nis{j};
end

% ----- Matched filters -----
ft = fliplr(p);
fs = fliplr(ps);

% Plot magnitude of filter responses
Fft = abs(fft(ft, 1024));
Ffs = abs(fft(fs, 1024));
figure;
semilogy(Fft); hold on; semilogy(Ffs);
legend('Square Pulse', 'Sinc Pulse');
title('Filter Spectral Responses (log scale)');

% ----- Filter the signals -----
e = cell(1,8); es = cell(1,8);
for j = 1:8
    e{j} = conv(d{j}, ft, 'same');
    es{j} = conv(ds{j}, fs, 'same');
end

% ----- Plot filtered signal (10 bits) -----
sample_points = (1:10) * bit_samples;
figure;
plot(e{1}); hold on;
stem(sample_points, e{1}(sample_points), 'r');
title('Filtered e^i(t) with Sample Points');

figure;
plot(es{1}); hold on;
stem(sample_points, es{1}(sample_points), 'g');
title('Filtered e^s(t) with Sample Points');

% ----- BER calculation -----
BER_i = zeros(1,8);
BER_s = zeros(1,8);

for j = 1:8
    ri = e{j}(bit_samples:bit_samples:(N*bit_samples));
    rj = es{j}(bit_samples:bit_samples:(N*bit_samples));
    decoded_i = ri > 0;
    decoded_s = rj > 0;

    BER_i(j) = sum(decoded_i ~= bb) / N;
    BER_s(j) = sum(decoded_s ~= bb) / N;
end

% ----- Plot BER vs. SNR -----
figure;
semilogy(SNR_dB, BER_i, '-o'); hold on;
semilogy(SNR_dB, BER_s, '-x');
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
legend('Square Pulse', 'Sinc Pulse');
title('BER vs. SNR');
grid on;
