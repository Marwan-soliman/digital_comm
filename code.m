%% ===================== COMMON SETUP =====================
clear; clc; close all;

N_sym   = 2e6;                   % symbols per branch / per scheme
E_bit   = 1;
EbN0_dB = 0:1:10;
EbN0    = 10.^(EbN0_dB/10);
N0      = E_bit ./ EbN0;
sigma   = sqrt(N0/2);            % noise std per real dimension

Q = @(x) 0.5*erfc(x/sqrt(2));

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
p = Q(sqrt(2*EbN0));
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