%% ===================== COMMON SETUP =====================
clear; clc; close all;

N_sym   = 2e6;                   % symbols per branch / per scheme
E_bit   = 1;
EbN0_dB = 0:1:10;
EbN0    = 10.^(EbN0_dB/10);
N0      = E_bit ./ EbN0;
sigma   = sqrt(N0/2);            % noise std per real dimension

%% ===================== PART 1: BPSK =====================
BER_bpsk = zeros(size(EbN0_dB));

bits_bpsk = randi([0 1], 1, N_sym);

for k = 1:length(EbN0_dB)
    % mapper
    s = (2*bits_bpsk - 1) * sqrt(E_bit);
    % channel
    r = s + sigma(k) * randn(size(s));
    % demapper
    bits_hat = (r > 0);
    % BER
    BER_bpsk(k) = sum(bits_hat ~= bits_bpsk) / N_sym;
    if EbN0_dB(k) == 10, r_bpsk10 = r; end                      % keep 10 dB received symbols (r is reused later)
end

%% ===================== PART 2: QPSK (Fig. 2, Gray) =====================
% Two independent BPSK branches: I and Q, each with its own bits and noise
BER_qpsk = zeros(size(EbN0_dB));

bits_I = randi([0 1], 1, N_sym);          % I-branch bits
bits_Q = randi([0 1], 1, N_sym);          % Q-branch bits

for k = 1:length(EbN0_dB)
    % ---- I branch ----
    sI   = (2*bits_I - 1) * sqrt(E_bit);              % mapper
    rI   = sI + sigma(k) * randn(size(sI));           % own noise
    hatI = (rI > 0);                                  % demapper

    % ---- Q branch ----
    sQ   = (2*bits_Q - 1) * sqrt(E_bit);              % mapper
    rQ   = sQ + sigma(k) * randn(size(sQ));           % own noise
    hatQ = (rQ > 0);                                  % demapper

    % BER over both bit streams (2*N_sym bits total)
    errs = sum(hatI ~= bits_I) + sum(hatQ ~= bits_Q);
    BER_qpsk(k) = errs / (2*N_sym);
    if EbN0_dB(k) == 10, tx_qpsk10 = sI + 1j*sQ; r_qpsk10 = rI + 1j*rQ; end   % sI/sQ/rI/rQ are reused in Part 3
end

%% ===================== PART 3: QPSK (Fig. 3, non-Gray) =====================
% Here the mapping itself couples the two bits, so the branches cannot be
% fully independent: I carries b1, Q carries (b1 XOR b2).
BER_qpsk_ng = zeros(size(EbN0_dB));

b1 = randi([0 1], 1, N_sym);
b2 = randi([0 1], 1, N_sym);

for k = 1:length(EbN0_dB)
    % mapper
    sI = (2*b1 - 1) * sqrt(E_bit);
    sQ = (2*xor(b1, b2) - 1) * sqrt(E_bit);

    % channel (independent noise on each dimension)
    rI = sI + sigma(k) * randn(size(sI));
    rQ = sQ + sigma(k) * randn(size(sQ));

    % demapper
    b1_hat = (rI > 0);
    x_hat  = (rQ > 0);                 % estimate of b1 XOR b2
    b2_hat = xor(b1_hat, x_hat);

    errs = sum(b1_hat ~= b1) + sum(b2_hat ~= b2);
    BER_qpsk_ng(k) = errs / (2*N_sym);
end

%% ===================== THEORY =====================
p = 0.5 * erfc(sqrt(EbN0));
BER_bpsk_th    = p;
BER_qpsk_th    = p;
BER_qpsk_ng_th = (3*p - 2*p.^2) / 2;

%% ===================== PLOT 1: BPSK vs QPSK =====================
figure;
semilogy(EbN0_dB, BER_bpsk, 'bo-', 'LineWidth', 1.2); hold on;
semilogy(EbN0_dB, BER_qpsk, 'rs-', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_bpsk_th, 'k--', 'LineWidth', 1.2);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BER vs E_b/N_0: BPSK and QPSK');
legend('BPSK (sim)', 'QPSK (sim)', 'Theory', 'Location', 'southwest');
ylim([1e-6 1]);

%% ===================== PLOT 2: QPSK Fig. 2 vs Fig. 3 =====================
figure;
semilogy(EbN0_dB, BER_qpsk,    'rs-', 'LineWidth', 1.2); hold on;
semilogy(EbN0_dB, BER_qpsk_ng, 'g^-', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_qpsk_th,    'k--', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_qpsk_ng_th, 'm--', 'LineWidth', 1.2);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('QPSK: Fig. 2 (Gray) vs Fig. 3 (non-Gray)');
legend('Fig. 2 (sim)', 'Fig. 3 (sim)', 'Fig. 2 (theory)', 'Fig. 3 (theory)', 'Location', 'southwest');
ylim([1e-6 1]);

%% ===================== PART 4: 8PSK (Fig. 2) =====================
BER_8psk = zeros(size(EbN0_dB));
N_sym_8psk = floor(N_sym / 3); 
bits_8psk = randi([0 1], 3, N_sym_8psk);

% Compact Gray Mapping Table (000, 001, 010, 011, 100, 101, 110, 111)
phase_map = [0, pi/4, 3*pi/4, pi/2, 7*pi/4, 3*pi/2, pi, 5*pi/4];
sym_idx_8psk = bits_8psk(1,:)*4 + bits_8psk(2,:)*2 + bits_8psk(3,:) + 1;

% Generate ideal constellation and transmit signal
const_8psk = sqrt(3 * E_bit) * exp(1j * phase_map); 
s_8psk = const_8psk(sym_idx_8psk);

for k = 1:length(EbN0_dB)
    % Channel
    r = s_8psk + sigma(k) * (randn(1, N_sym_8psk) + 1j * randn(1, N_sym_8psk));
    
    % Demapper (Minimum Euclidean Distance)
    [~, est_idx] = min(abs(r.' - const_8psk), [], 2);
    est_bits = dec2bin(est_idx - 1, 3) - '0';
    
    % BER calculation
    BER_8psk(k) = sum(sum(bits_8psk' ~= est_bits)) / (3 * N_sym_8psk);
    if EbN0_dB(k) == 10, tx_8psk10 = s_8psk; rx_8psk10 = r; end
end

%% ===================== PART 5: 16-QAM (Fig. 2) =====================
BER_16qam = zeros(size(EbN0_dB));
N_sym_16qam = floor(N_sym / 4);

bits_I_16qam = randi([0 1], 2, N_sym_16qam); 
bits_Q_16qam = randi([0 1], 2, N_sym_16qam); 

% Compact 4-PAM Mapper (00->-3, 01->-1, 10->+3, 11->+1)
pam_map = [-3, -1, 3, 1];
idx_I = bits_I_16qam(1,:)*2 + bits_I_16qam(2,:) + 1;
idx_Q = bits_Q_16qam(1,:)*2 + bits_Q_16qam(2,:) + 1;

a = sqrt(4 * E_bit / 10); % Scale factor for Es = 4*Eb
s_I = a * pam_map(idx_I);
s_Q = a * pam_map(idx_Q);

for k = 1:length(EbN0_dB)
    % Channel
    r_I = s_I + sigma(k) * randn(size(s_I));
    r_Q = s_Q + sigma(k) * randn(size(s_Q));
    
    % Demapper
    pam_demap = @(r) [r > 0; abs(r/a) < 2]; 
    hat_I = pam_demap(r_I);
    hat_Q = pam_demap(r_Q);
    
    % BER calculation
    errs = sum(sum(hat_I ~= bits_I_16qam)) + sum(sum(hat_Q ~= bits_Q_16qam));
    BER_16qam(k) = errs / (4 * N_sym_16qam);
    if EbN0_dB(k) == 10, tx_16qam10 = s_I + 1j*s_Q; rx_16qam10 = r_I + 1j*r_Q; end
end

%% ===================== THEORY (8PSK & 16QAM) =====================
BER_8psk_th  = (1/3) * erfc(sqrt(3 * EbN0) * sin(pi/8));
BER_16qam_th = (3/8) * erfc(sqrt(0.4 * EbN0));

%% ===================== PART 6: Coherent BFSK =====================
% 2D orthogonal signal space:
%   bit 0 -> s1 = [sqrt(E_bit); 0] along phi_1
%   bit 1 -> s2 = [0; sqrt(E_bit)] along phi_2
BER_bfsk  = zeros(size(EbN0_dB));
bits_bfsk = randi([0 1], 1, N_sym);

for k = 1:length(EbN0_dB)
    % mapper (projections onto phi_1 and phi_2)
    s1 = (bits_bfsk == 0) * sqrt(E_bit);
    s2 = (bits_bfsk == 1) * sqrt(E_bit);
    
    % channel (independent AWGN with variance N0/2 per orthogonal dimension)
    r1 = s1 + sigma(k) * randn(size(s1));
    r2 = s2 + sigma(k) * randn(size(s2));
    
    % demapper (ML decision rule: decide 1 if r2 > r1, else 0)
    bits_hat = (r2 > r1);
    
    % BER calculation
    BER_bfsk(k) = sum(bits_hat ~= bits_bfsk) / N_sym;
    if EbN0_dB(k) == 10, tx_bfsk10 = s1 + 1j*s2; rx_bfsk10 = r1 + 1j*r2; end
end

%% ===================== THEORY (BFSK) =====================
BER_bfsk_th  = 0.5 * erfc(sqrt(0.5 * EbN0));

%% ===================== PART 7: BFSK Baseband PSD Simulation =====================
Tb    = 1;                                        % normalized bit duration (s)
Ns    = 64;                                       % oversampling factor (samples/bit)
Fs    = Ns / Tb;                                  % sampling frequency (Hz)                                    

% Time vector for one bit interval [0, Tb)
t_bit = (0:Ns-1) / Fs;

% Baseband equivalent signals (Lecture convention: fc = f1, Delta_f = 1/Tb)
s1_bb = sqrt(2*E_bit/Tb) * ones(1, Ns);
s2_bb = sqrt(2*E_bit/Tb) * exp(1j * 2 * pi * (1/Tb) * t_bit);

% Generate random bits and synthesize continuous-time baseband stream
bits_psd = randi([0 1], N_sym, 1);
tx_mat   = (bits_psd == 0) * s1_bb + (bits_psd == 1) * s2_bb;
x_bb     = reshape(tx_mat.', 1, []);

% Estimate centered two-sided PSD via Welch's method
win_len = 4096;
overlap = win_len / 2;
nfft    = 8192;
[S_B, f] = pwelch(x_bb, hamming(win_len), overlap, nfft, Fs, 'centered');

%% ===================== PLOT 4: Baseband PSD =====================

% Linear Normalized Baseband PSD
figure;
plot(f * Tb, S_B / E_bit, 'b-', 'LineWidth', 1.4);
grid on;
xlabel('Normalized Frequency (f \cdot T_b)');
ylabel('Normalized PSD, S_B(f) / E_b');
title('Simulated Baseband Power Spectral Density of Coherent BFSK');
xlim([-3 4]);

% dB Normalized Baseband PSD
figure;
plot(f * Tb, 10*log10(S_B / E_bit), 'b-', 'LineWidth', 1.2);
grid on;
xlabel('Normalized Frequency (f \cdot T_b)');
ylabel('Normalized PSD, 10 log_{10}(S_B(f) / E_b) (dB)');
title('Simulated Baseband PSD of Coherent BFSK (dB Scale)');
xlim([-3 4]); ylim([-50 5]);

%% ===================== PLOT 5: Integrated Comparison =====================
figure;
% Simulation curves
semilogy(EbN0_dB, BER_bpsk,  'bo-', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
semilogy(EbN0_dB, BER_qpsk,  'rs-', 'LineWidth', 1.2, 'MarkerSize', 6);
semilogy(EbN0_dB, BER_8psk,  'g^-', 'LineWidth', 1.2, 'MarkerSize', 6);
semilogy(EbN0_dB, BER_16qam, 'md-', 'LineWidth', 1.2, 'MarkerSize', 6);
semilogy(EbN0_dB, BER_bfsk,  'kv-', 'LineWidth', 1.2, 'MarkerSize', 6);

% Theoretical curves
semilogy(EbN0_dB, BER_bpsk_th,  'b--', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_qpsk_th,  'r--', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_8psk_th,  'g--', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_16qam_th, 'm--', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_bfsk_th,  'k--', 'LineWidth', 1.2);

grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('Integrated BER vs E_b/N_0 for BPSK, QPSK, 8PSK, 16QAM, and BFSK');
legend('BPSK (sim)', 'QPSK (sim)', '8PSK (sim)', '16QAM (sim)', 'BFSK (sim)', ...
       'BPSK (th)', 'QPSK (th)', '8PSK (th)', '16QAM (th)', 'BFSK (th)', ...
       'Location', 'southwest');
ylim([1e-6 1]);

%% ===================== EXTRA: Constellations @ 10 dB =====================
n   = 1:5000;                                                % Extract first 5000 symbols at 10 dB
TX  = {s, tx_qpsk10, tx_8psk10, tx_16qam10, tx_bfsk10};      % Transmitted symbols array
RX  = {r_bpsk10, r_qpsk10, rx_8psk10, rx_16qam10, rx_bfsk10};% Received symbols array
ttl = {'BPSK', 'QPSK', '8PSK', '16QAM', 'BFSK'};
figure;
for m = 1:5
    subplot(2,3,m);
    plot(real(RX{m}(n)), imag(RX{m}(n)), '.', 'MarkerSize', 4); hold on;   % Received points
    plot(real(TX{m}(n)), imag(TX{m}(n)), 'r+', 'MarkerSize', 8);           % Ideal targets
    grid on; axis equal; title([ttl{m} ' @ E_b/N_0 = 10 dB']);
    if m == 5
        xlabel('\phi_1'); ylabel('\phi_2');
    else
        xlabel('I'); ylabel('Q');
    end
end
