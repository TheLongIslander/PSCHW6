RUID = 208001821;
rng(RUID);

T = 2;                 
A = 1;                 
Ts = 0.02;           
R = 1 / T;
N = 100000;
bit_samples = round(T / Ts);
L = N * bit_samples;

t_sinc = -5*T:Ts:5*T;
ps = A * sinc(t_sinc / T);
ps_delay = floor(length(ps)/2);

% Bit sequence
bb = randi([0, 1], 1, N);
disp('First 10 bits:');
disp(bb(1:10));

% Define square pulse
p = A * ones(1, bit_samples);
Ep = sum(p.^2) * Ts;
Eps = sum(ps.^2) * Ts;

% Preallocate signals
s = zeros(1, L + length(ps));            % square pulse
ss = zeros(1, L + 2*ps_delay);           % sinc pulse
valid = false(1, N);                     % track valid bits for sinc

% Modulate both signals
for i = 1:N
    idx = (i-1)*bit_samples + 1;
    center_idx = idx + bit_samples/2 - ps_delay;
    polarity = 2 * bb(i) - 1;

    s(idx:idx+bit_samples-1) = s(idx:idx+bit_samples-1) + polarity * p;

    if center_idx > 0 && (center_idx + length(ps) - 1) <= length(ss)
        ss(center_idx:center_idx+length(ps)-1) = ...
            ss(center_idx:center_idx+length(ps)-1) + polarity * ps;
        valid(i) = true;
    end
end

% Plot pulses
figure;
subplot(2,1,1); plot(p); title('Square Pulse p(t)');
subplot(2,1,2); plot(ps); title('Sinc Pulse ps(t)');

fprintf('Energy of p(t): %.4f, Energy of ps(t): %.4f\n', Ep, Eps);
fprintf('Sinc modulated %d out of %d bits (%.2f%%)\n', sum(valid), N, 100*sum(valid)/N);

% Noise generation
SNR_dB = 0:1:7;
SNR = 10.^(SNR_dB / 10);
ni = cell(1, 8); nis = cell(1, 8);
for j = 1:8
    vi = Ep / SNR(j);
    vis = Eps / SNR(j);
    ni{j} = sqrt(vi) * randn(1, length(s));
    nis{j} = sqrt(vis) * randn(1, length(ss));
end

% Add noise
d = cell(1, 8); ds = cell(1, 8);
for j = 1:8
    d{j} = s + ni{j};
    ds{j} = ss + nis{j};
end

% Matched filters
ft = fliplr(p);
fs = fliplr(ps);

% Frequency response
Fft = abs(fft(ft, 1024));
Ffs = abs(fft(fs, 1024));
figure;
semilogy(Fft); hold on; semilogy(Ffs);
legend('Square Pulse', 'Sinc Pulse');
title('Filter Spectral Responses (log scale)');

% Filtering
e = cell(1,8); es = cell(1,8);
for j = 1:8
    e{j} = conv(d{j}, ft, 'same');
    es{j} = conv(ds{j}, fs, 'same');
end

% Sample points (bit centers)
sample_points_i = bit_samples/2 : bit_samples : N * bit_samples - bit_samples/2;
sample_points_s = sample_points_i - ps_delay;
sample_points_s = sample_points_s(sample_points_s > 0);  % ✅ Filter invalid indices

% Sanity check plots
figure;
plot(e{1}); hold on;
stem(sample_points_i(1:10), e{1}(sample_points_i(1:10)), 'r');
title('Filtered e^i(t) with Sample Points');

figure;
plot(es{1}); hold on;
stem(sample_points_s(1:10), es{1}(sample_points_s(1:10)), 'g');
title('Filtered e^s(t) with Sample Points');

% BER Calculation
BER_i = zeros(1,8);
BER_s = zeros(1,8);

for j = 1:8
    ri = e{j}(sample_points_i);
    rj = es{j}(sample_points_s);
    decoded_i = ri > 0;
    decoded_s = rj > 0;

    % Truncate bb to match number of decoded bits
    BER_i(j) = sum(decoded_i ~= bb(1:length(decoded_i))) / length(decoded_i);
    BER_s(j) = sum(decoded_s ~= bb(1:length(decoded_s))) / length(decoded_s);
end

% BER plot
figure;
semilogy(SNR_dB, BER_i, '-o'); hold on;
semilogy(SNR_dB, BER_s, '-x');
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
legend('Square Pulse', 'Sinc Pulse');
title('BER vs. SNR');
grid on;
