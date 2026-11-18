%% ===================== COMMON SETUP =====================
clear; clc; close all;

N_sym   = 2e6;                   % symbols per SNR point
E_bit   = 1;
EbN0_dB = 0:2:26;
EbN0    = 10.^(EbN0_dB/10);
N0      = E_bit ./ EbN0;
sigma   = sqrt(N0/2);            % noise std per real dimension

%% ===================== PART 1: BPSK (Rayleigh Flat Fading) =====================
BER_bpsk_ray = zeros(size(EbN0_dB));
bits_bpsk    = randi([0 1], 1, N_sym);

% Mapper (0 -> -sqrt(Eb), 1 -> +sqrt(Eb))
x_bpsk = (2*bits_bpsk - 1) * sqrt(E_bit);

for k = 1:length(EbN0_dB)
    % Rayleigh flat fading channel & AWGN
    h = sqrt(0.5) * (randn(1, N_sym) + 1j*randn(1, N_sym));
    n = sigma(k)  * (randn(1, N_sym) + 1j*randn(1, N_sym));

    % Received signal
    y = h .* x_bpsk + n;

    % Channel inversion & hard decision
    y_eq     = y ./ h;
    bits_hat = (real(y_eq) > 0);

    % BER calculation
    BER_bpsk_ray(k) = sum(bits_hat ~= bits_bpsk) / N_sym;

    if EbN0_dB(k) == 20
        tx_bpsk20 = x_bpsk;
        rx_bpsk20 = y_eq;
    end
end

%% ===================== PART 2: QPSK (Rayleigh Flat Fading) =====================
BER_qpsk_ray = zeros(size(EbN0_dB));
bits_I       = randi([0 1], 1, N_sym);
bits_Q       = randi([0 1], 1, N_sym);

% Mapper (Gray: 00->-1-j, 01->-1+j, 10->+1-j, 11->+1+j)
sI     = (2*bits_I - 1) * sqrt(E_bit);
sQ     = (2*bits_Q - 1) * sqrt(E_bit);
x_qpsk = sI + 1j*sQ;

for k = 1:length(EbN0_dB)
    % Rayleigh flat fading channel & AWGN
    h = sqrt(0.5) * (randn(1, N_sym) + 1j*randn(1, N_sym));
    n = sigma(k)  * (randn(1, N_sym) + 1j*randn(1, N_sym));

    % Received signal
    y = h .* x_qpsk + n;

    % Channel inversion & hard decision
    y_eq = y ./ h;
    hatI = (real(y_eq) > 0);
    hatQ = (imag(y_eq) > 0);

    % BER calculation
    errs = sum(hatI ~= bits_I) + sum(hatQ ~= bits_Q);
    BER_qpsk_ray(k) = errs / (2*N_sym);

    if EbN0_dB(k) == 20
        tx_qpsk20 = x_qpsk;
        rx_qpsk20 = y_eq;
    end
end

%% ===================== PART 3: 16-QAM (Rayleigh Flat Fading) =====================
BER_16qam_ray = zeros(size(EbN0_dB));
bits_I_16qam  = randi([0 1], 2, N_sym);
bits_Q_16qam  = randi([0 1], 2, N_sym);

% Mapper (Gray 4-PAM per axis: 00->-3, 01->-1, 10->+3, 11->+1)
pam_map = [-3, -1, 3, 1];
idx_I   = bits_I_16qam(1,:)*2 + bits_I_16qam(2,:) + 1;
idx_Q   = bits_Q_16qam(1,:)*2 + bits_Q_16qam(2,:) + 1;

a       = sqrt(4 * E_bit / 10);  % scale factor for Es = 4*Eb
s_I     = a * pam_map(idx_I);
s_Q     = a * pam_map(idx_Q);
x_16qam = s_I + 1j*s_Q;

pam_demap = @(r) [r > 0; abs(r/a) < 2];

for k = 1:length(EbN0_dB)
    % Rayleigh flat fading channel & AWGN
    h = sqrt(0.5) * (randn(1, N_sym) + 1j*randn(1, N_sym));
    n = sigma(k)  * (randn(1, N_sym) + 1j*randn(1, N_sym));

    % Received signal
    y = h .* x_16qam + n;

    % Channel inversion & hard decision
    y_eq  = y ./ h;
    hat_I = pam_demap(real(y_eq));
    hat_Q = pam_demap(imag(y_eq));

    % BER calculation
    errs = sum(sum(hat_I ~= bits_I_16qam)) + sum(sum(hat_Q ~= bits_Q_16qam));
    BER_16qam_ray(k) = errs / (4*N_sym);

    if EbN0_dB(k) == 20
        tx_16qam20 = x_16qam;
        rx_16qam20 = y_eq;
    end
end

%% ===================== THEORY =====================
f = @(c) 0.5 * (1 - sqrt((c * EbN0) ./ (1 + c * EbN0)));

BER_bpsk_ray_th  = f(1);
BER_qpsk_ray_th  = f(1);
BER_16qam_ray_th = (3/4)*f(0.4) + (1/2)*f(3.6) - (1/4)*f(10);

% AWGN theoretical references
BER_bpsk_awgn_th  = 0.5 * erfc(sqrt(EbN0));
BER_qpsk_awgn_th  = 0.5 * erfc(sqrt(EbN0));
BER_16qam_awgn_th = (3/8) * erfc(sqrt(0.4 * EbN0));

%% ===================== PLOT 1: Integrated Rayleigh BER Comparison =====================
figure;
semilogy(EbN0_dB, BER_bpsk_ray,  'bo-', 'LineWidth', 1.2, 'MarkerSize', 6); hold on;
semilogy(EbN0_dB, BER_qpsk_ray,  'rs-', 'LineWidth', 1.2, 'MarkerSize', 6);
semilogy(EbN0_dB, BER_16qam_ray, 'md-', 'LineWidth', 1.2, 'MarkerSize', 6);

semilogy(EbN0_dB, BER_bpsk_ray_th,  'b--', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_qpsk_ray_th,  'r--', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_16qam_ray_th, 'm--', 'LineWidth', 1.2);

grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BER vs E_b/N_0 over Rayleigh Flat Fading Channel');
legend('BPSK (sim)', 'QPSK (sim)', '16-QAM (sim)', ...
       'BPSK (th)', 'QPSK (th)', '16-QAM (th)', 'Location', 'southwest');
ylim([1e-4 1]);

%% ===================== EXTRA - PLOT 2: Rayleigh Fading vs AWGN =====================
figure;
semilogy(EbN0_dB, BER_bpsk_ray,     'bo-', 'LineWidth', 1.2, 'MarkerSize', 5); hold on;
semilogy(EbN0_dB, BER_qpsk_ray,     'rs-', 'LineWidth', 1.2, 'MarkerSize', 5);
semilogy(EbN0_dB, BER_16qam_ray,    'md-', 'LineWidth', 1.2, 'MarkerSize', 5);
semilogy(EbN0_dB, BER_bpsk_awgn_th, 'k-.', 'LineWidth', 1.2);
semilogy(EbN0_dB, BER_16qam_awgn_th,'g-.', 'LineWidth', 1.2);

grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title('BER Comparison: Rayleigh Flat Fading vs AWGN');
legend('BPSK Rayleigh (sim)', 'QPSK Rayleigh (sim)', '16-QAM Rayleigh (sim)', ...
       'BPSK/QPSK AWGN (th)', '16-QAM AWGN (th)', 'Location', 'southwest');
ylim([1e-5 1]);

%% ===================== EXTRA - PLOT 3: Equalized Constellations @ 20 dB =====================
n_plot = 1:5000;
TX     = {tx_bpsk20, tx_qpsk20, tx_16qam20};
RX     = {rx_bpsk20, rx_qpsk20, rx_16qam20};
ttl    = {'BPSK', 'QPSK', '16-QAM'};

figure;
for m = 1:3
    subplot(1,3,m);
    plot(real(RX{m}(n_plot)), imag(RX{m}(n_plot)), '.', 'MarkerSize', 4); hold on;
    plot(real(TX{m}(n_plot)), imag(TX{m}(n_plot)), 'r+', 'MarkerSize', 8, 'LineWidth', 1.2);
    grid on; axis equal; xlim([-4 4]); ylim([-4 4]);
    xlabel('I'); ylabel('Q');
    title([ttl{m} ' Equalized @ 20 dB']);
end