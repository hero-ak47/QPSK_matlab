clc; clear; close all

%% SECTION 1: BẮT ĐẦU THU TÍN HIỆU
Fs = 48000;
CAPTURE_SEC = 7;

disp('Đang ghi âm 7 giây...');
r = audiorecorder(Fs, 24, 1);

recordblocking(r, CAPTURE_SEC);

disp('Đã dừng ghi âm. Đang xử lý tín hiệu...');
y = getaudiodata(r, 'double');
% test
% ================= NẠP THÔNG TIN MÃ HÓA TỪ TX (MỚI THÊM) =================
load('tx_coding_info.mat');   % nap: trellis, n_info, n_tail, info_bits,
                               %      constraint_length, code_generators,
                               %      interleaver_pattern
fprintf('=== DA NAP THONG TIN MA HOA TU TX ===\n');
fprintf('Constraint length K = %d, generator = [%d %d] (octal)\n', ...
        constraint_length, code_generators(1), code_generators(2));
fprintf('So bit thong tin goc: %d, so bit tail: %d\n\n', n_info, n_tail);

% 1. THÔNG SỐ VÀ TÁI TẠO PREAMBLE LÀM CHUẨN MẪU
Fs = 48000; Fc = 12000;
Rs = 200;   % symbol rate
L = Fs/Rs;  % samples/symbol

% RRC
beta = 0.5; span = 6;
h_rrc = rcosdesign(beta, span, L, 'sqrt');

% Zadoff-Chu
Nzc = 63;
u   = 25;
n   = (0:Nzc-1).';
pss = exp(-1j*pi*u*n.*(n+1)/Nzc);

% preamble la ZC de dong bo
preamble_sym_rx = pss;
preamble_bb_ideal = upfirdn(preamble_sym_rx, h_rrc, L);

% ================= NẠP TÍN HIỆU =================
% [y, Fs_read] = audioread('tx_signal_1.wav');
% if Fs_read ~= Fs
%     warning('Fs trong file (%d) khác Fs khai báo (%d)!', Fs_read, Fs);
% end
% disp('Đã nạp tín hiệu từ tx_signal_1.wav');

% Khu offset
y = y - mean(y);

rayleighChan = comm.RayleighChannel( ...
    'SampleRate', 48000, ...
    'PathDelays', [0 0.0001 0.0002], ...
    'AveragePathGains', [0 -9 -10], ...
    'MaximumDopplerShift', 0);
y = rayleighChan(y);

% mo phong kenh truyen : CFO
CFO_true = 1;
y = y .* exp(1j*2*pi*CFO_true*(0:length(y)-1)'/Fs);
%awgn
SNR_dB = 5;
y = awgn(y, SNR_dB);
plot(y);

% 2. TIỀN XỬ LÝ: LỌC BĂNG THÔNG & HẠ TẦN XUỐNG BASEBAND
[b, a] = butter(5, [(Fc - 1000) (Fc + 1000)]/(Fs/2), 'bandpass');
y_filt = filter(b, a, y);

t_rx = (0:length(y_filt)-1)' / Fs;
y_bb = y_filt .* exp(-1j * 2 * pi * Fc * t_rx);

y_mf = conv(y_bb, h_rrc, 'same');
% ================= XUAT CAPTURE.BIN CHO C =================

% Dam bao dung 7 giay
N_expected = Fs * CAPTURE_SEC;

if length(y_mf) < N_expected
    y_mf(end+1:N_expected) = 0;
elseif length(y_mf) > N_expected
    y_mf = y_mf(1:N_expected);
end

% Tach I/Q
I = real(y_mf);
Q = imag(y_mf);

% Normalize chung theo peak
scale = max(abs([I; Q]));

if scale == 0
    scale = 1;
end

I_norm = I / scale;
Q_norm = Q / scale;

% Q1.31
I_q31 = int32(round(I_norm * 2147483647));
Q_q31 = int32(round(Q_norm * 2147483647));

% Interleave: I0 Q0 I1 Q1 ...
iq = zeros(2*N_expected, 1, 'int32');
iq(1:2:end) = I_q31;
iq(2:2:end) = Q_q31;

fid = fopen('capture_1.bin', 'wb');
fwrite(fid, iq, 'int32');
fclose(fid);

fprintf('\n=== XUAT CAPTURE.BIN ===\n');
fprintf('Fs              = %d Hz\n', Fs);
fprintf('Thoi gian       = %d s\n', CAPTURE_SEC);
fprintf('So complex samp = %d\n', N_expected);
fprintf('So int32        = %d\n', length(iq));
fprintf('File size       = %.2f MB\n', ...
    dir('capture_1.bin').bytes / 1024^2);
%-----------------------------------------------------
% 3. ĐỒNG BỘ KHUNG (FRAME SYNCHRONIZATION)
Nzc       = 63;
N_guard   = 14;
N_payload = 500;
N_between_symbols = Nzc + N_guard + N_payload;
T0 = N_between_symbols * L / Fs;

[xc, lags] = xcorr(y_mf, preamble_bb_ideal);
xc_abs = abs(xc);
noise_floor = median(xc_abs);

search_limit_lag = round(length(y_mf) / 2);
valid_idx1 = find(lags >= 0 & lags <= search_limit_lag);
if isempty(valid_idx1)
    error('Vùng tìm kiếm đỉnh 1 không hợp lệ -> kiểm tra lại độ dài tín hiệu!');
end

[peak1_val, rel_idx1] = max(xc_abs(valid_idx1));
idx1 = valid_idx1(rel_idx1);
lag1 = lags(idx1);
start_sample = lag1 + 1;

fprintf('\n--- DEBUG: FRAME SYNC (đỉnh 1) ---\n');
fprintf('start_sample = %d (length y_mf = %d)\n', start_sample, length(y_mf));
fprintf('Peak1/Noise = %.2f\n', peak1_val/noise_floor);
if peak1_val/noise_floor < 5
    warning('Peak1/Noise thấp -> có thể sync SAI vị trí!');
end

if start_sample < 1 || (start_sample + (Nzc+N_guard+N_payload)*L) > length(y_mf)
    error('Không tìm thấy khung tín hiệu nguyên vẹn. Hãy thử thu phát lại!');
end

%------------------------------------------------------------------
% ƯỚC LƯỢNG DOPPLER: TÌM TIẾP ĐỈNH 2 (preamble cuối) — CHỈ trong nửa sau
% (MỚI THÊM LẠI: TX bây giờ đã có preamble thứ 2 ở cuối frame: [...; pss])
%------------------------------------------------------------------
fprintf('\n--- DEBUG: 2-PILOT DOPPLER ESTIMATION ---\n');
T0 = N_between_symbols * L / Fs;
fprintf('T0 (khoảng cách danh định giữa 2 preamble) = %.6f s\n', T0);

valid_idx2 = find(lags > search_limit_lag & lags <= length(y_mf)-1);
if isempty(valid_idx2)
    error('Vùng tìm kiếm đỉnh 2 không hợp lệ -> kiểm tra lại độ dài tín hiệu!');
end

[peak2_val, rel_idx2] = max(xc_abs(valid_idx2));
idx2 = valid_idx2(rel_idx2);
lag2 = lags(idx2);

fprintf('Đỉnh 2 (preamble CUỐI, nửa sau tín hiệu): lag = %d mẫu, biên độ = %.3f\n', lag2, peak2_val);
fprintf('Peak2/Noise = %.2f\n', peak2_val/noise_floor);
if peak2_val/noise_floor < 5
    warning('Đỉnh 2 yếu -> có thể xác định SAI vị trí preamble cuối!');
end

% ---- Tính hệ số Doppler (time-scale) ----
sample_diff  = lag2 - lag1;
T_rx         = sample_diff / Fs;
a_hat_2pilot = T0/T_rx - 1;

fprintf('T_rx (đo được)   = %.6f s\n', T_rx);
fprintf('a_hat (2-pilot)  = %.6f\n', a_hat_2pilot);
fprintf('f_D tương đương ~ %.3f Hz (tại Fc = %d Hz)\n', a_hat_2pilot*Fc, Fc);

% ---- ÁP DỤNG BÙ DOPPLER TIME-SCALE ----
y = doppler_resample(y - mean(y), a_hat_2pilot, Fs);
fprintf('Đã bù Doppler time-scale, độ dài tín hiệu sau bù: %d mẫu\n', length(y));

% ---- LẶP LẠI TIỀN XỬ LÝ + ĐỒNG BỘ KHUNG TRÊN TÍN HIỆU ĐÃ BÙ ----
y_filt = filter(b, a, y);
t_rx   = (0:length(y_filt)-1)' / Fs;
y_bb   = y_filt .* exp(-1j * 2 * pi * Fc * t_rx);
y_mf   = conv(y_bb, h_rrc, 'same');

[xc, lags] = xcorr(y_mf, preamble_bb_ideal);
xc_abs = abs(xc);
noise_floor = median(xc_abs);

search_limit_lag = round(length(y_mf) / 2);
valid_idx1 = find(lags >= 0 & lags <= search_limit_lag);
[peak1_val, rel_idx1] = max(xc_abs(valid_idx1));
idx1 = valid_idx1(rel_idx1);
start_sample = lags(idx1) + 1;

fprintf('\n--- FRAME SYNC LẠI sau khi bù Doppler ---\n');
fprintf('start_sample (mới) = %d\n', start_sample);
fprintf('Peak1/Noise = %.2f\n', peak1_val/noise_floor);

if start_sample < 1 || (start_sample + (Nzc+N_guard+N_payload)*L) > length(y_mf)
    error('Không tìm thấy khung tín hiệu nguyên vẹn sau khi bù Doppler!');
end

% 4. TRÍCH XUẤT SYMBOL (DOWNSAMPLING)
total_sym_len = 577 + 63;
skip = 77;
rx_symbols = zeros(total_sym_len, 1);

first_sym_offset = (span*L/2) + 1;

for k = 1:total_sym_len
    idx = start_sample + first_sym_offset + (k-1)*L - 1;
    rx_symbols(k) = y_mf(idx);
end

rx_payload = rx_symbols(skip+1:end - Nzc);

% ---- Vị trí pilot & chuỗi pilot ZC ----
pilot_indices = 1:5:500;
data_indices  = setdiff(1:500, pilot_indices);

N_pilots = length(pilot_indices);
u_p = 7;
n_p = (0:N_pilots-1).';
zc_pilots = sqrt(2) * exp(-1j * pi * u_p * n_p .* (n_p + 1) / N_pilots);

%-------------------------------------------------------------------------
% 5. ƯỚC LƯỢNG & BÙ CFO BẰNG LEAST-SQUARES TRÊN PHA PILOT ZC
%-------------------------------------------------------------------------
[cfo_est, rx_payload_cfo, info] = estimate_and_correct_cfo(...
    rx_payload, pilot_indices, zc_pilots, L, Fs, ...
    'CFO_true', CFO_true, 'Verbose', true, 'PlotConst', true);

%----------------------------------------------------------------
% 6. CÂN BẰNG KÊNH BẰNG RLS (DÙNG PILOT ZC)
%----------------------------------------------------------------
[rx_payload_eq, rx_data_eq, W_final, W_track, info] = rls_equalize_zc(...
    rx_payload_cfo, pilot_indices, zc_pilots, data_indices, ...
    'PlotConst', true, 'Verbose', true);

% ---- Bù lệch pha dư (residual phase offset) ----
constellation = [pi/4, 3*pi/4, -3*pi/4, -pi/4];
rx_phase = angle(rx_data_eq(1:400));
delta_rx = zeros(1,400);
for i = 1:400
    phase_error = angle(exp(1j*(rx_phase(i) - constellation)));
    [~, idx] = min(abs(phase_error));
    delta_rx(i) = phase_error(idx);
end
phi = angle(mean(exp(1j*delta_rx)));
rx_data_eq = rx_data_eq .* exp(-1j*phi);

fprintf('\nResidual phase offset sau RLS = %.4f rad (%.2f do)\n', phi, rad2deg(phi));

% 7. GIẢI ĐIỀU CHẾ QPSK (DEMAPPING) -- RA BIT DA MA HOA VA DA INTERLEAVE
rx_coded_bits_il = zeros(length(rx_data_eq)*2, 1);
rx_coded_bits_il(1:2:end) = real(rx_data_eq) < 0;
rx_coded_bits_il(2:2:end) = imag(rx_data_eq) < 0;

% ================= DEINTERLEAVER (MỚI THÊM) =================
% Xao tron NGUOC LAI theo dung interleaver_pattern da nap tu tx_coding_info.mat
rx_coded_bits = zeros(size(rx_coded_bits_il));
rx_coded_bits(interleaver_pattern) = rx_coded_bits_il;

fprintf('Da ap dung de-interleaver (khoi phuc thu tu bit goc)\n');

% ================= 8. GIẢI MÃ KÊNH: VITERBI DECODE (MỚI THÊM) =================
fprintf('\n=== GIAI MA KENH (VITERBI) ===\n');

% ---- BER TRUOC KHI GIAI MA (tren bit da encode, de so sanh) ----
coded_bits_ideal = convenc([info_bits; zeros(n_tail,1)], trellis);
num_errors_raw = sum(rx_coded_bits ~= coded_bits_ideal);
BER_before_decode = num_errors_raw / length(coded_bits_ideal);

fprintf('BER TRUOC khi decode (tren bit da ma hoa, kenh loi truc tiep): %f\n', BER_before_decode);

% ---- Viterbi decode ----
tblen = 5 * constraint_length;      % traceback depth khuyen nghi ~5*K = 35
decoded_bits = vitdec(rx_coded_bits, trellis, tblen, 'term', 'hard');

% ---- Bo 6 bit tail, lay lai bit thong tin goc ----
decoded_info = decoded_bits(1:n_info);

num_errors_after = sum(decoded_info ~= info_bits);
BER_after_decode = num_errors_after / n_info;

fprintf('BER SAU khi decode (tren bit thong tin goc, da sua loi):     %f\n', BER_after_decode);
fprintf('So bit loi: TRUOC = %d / %d  ->  SAU = %d / %d\n', ...
        num_errors_raw, length(coded_bits_ideal), num_errors_after, n_info);

if BER_after_decode < BER_before_decode
    fprintf('>> Coding CO HIEU QUA: giam BER tu %.4f xuong %.4f\n', ...
            BER_before_decode, BER_after_decode);
elseif BER_after_decode == 0
    fprintf('>> Coding SUA LOI HOAN TOAN (BER sau decode = 0)\n');
else
    fprintf('>> CANH BAO: Coding KHONG cai thien (co the SNR qua thap, vuot kha nang sua loi cua code)\n');
end

% 9. HIỂN THỊ KẾT QUẢ TỔNG QUAN
disp('--- KẾT QUẢ TRUYỀN NHẬN (CFO + RLS + Convolutional/Viterbi) ---');
fprintf('So bit thong tin goc: %d\n', n_info);
fprintf('So bit da ma hoa (truyen di): %d\n', length(coded_bits_ideal));
fprintf('BER truoc decode (raw channel): %f\n', BER_before_decode);
fprintf('BER sau decode (final):         %f\n', BER_after_decode);

figure;
scatter(real(rx_data_eq), imag(rx_data_eq), 'b.'); hold on;
scatter([1 -1 1 -1], [1 1 -1 -1], 'rx', 'LineWidth', 2);
title('Chòm sao tín hiệu nhận được (CFO + RLS, trước Viterbi decode)');
xlabel('In-phase'); ylabel('Quadrature');
grid on; axis square;

function y_out = doppler_resample(y, a_try, Fs, fc)
    if nargin < 4 || isempty(fc)
        fc = 0;
    end
    N = length(y);
    t_orig = (0:N-1)' / Fs;
    t_new = t_orig / (1 + a_try);
    y_resampled = interp1(t_orig, y, t_new, 'spline', 'extrap');
    if fc ~= 0
        y_out = y_resampled .* exp(-1i * 2 * pi * fc * a_try * t_orig);
    else
        y_out = y_resampled;
    end
end

function [cfo_est, rx_payload_cfo, info] = estimate_and_correct_cfo(rx_payload, pilot_indices, zc_pilots, L, Fs, varargin)
    p = inputParser;
    addParameter(p, 'CFO_true', []);
    addParameter(p, 'Verbose', true);
    addParameter(p, 'PlotConst', false);
    parse(p, varargin{:});
    CFO_true  = p.Results.CFO_true;
    verbose   = p.Results.Verbose;
    do_plot   = p.Results.PlotConst;

    rx_payload    = rx_payload(:);
    pilot_indices = pilot_indices(:);
    zc_pilots     = zc_pilots(:);

    rx_pilots_raw = rx_payload(pilot_indices);
    pilot_phase_raw = angle(rx_pilots_raw ./ zc_pilots);

    if verbose
        fprintf('\n--- DEBUG: PILOT PHASE (TRUOC UNWRAP) ---\n');
        disp(pilot_phase_raw(1:min(10, end)).');
    end

    pilot_phase_unwrapped = unwrap(pilot_phase_raw);

    idx_col = pilot_indices - 1;
    A = [ones(length(idx_col), 1), idx_col];

    coeffs = A \ pilot_phase_unwrapped;
    slope  = coeffs(2);

    Tsym    = L / Fs;
    cfo_est = slope / (2*pi*Tsym);

    phase_fit     = A * coeffs;
    residual      = pilot_phase_unwrapped - phase_fit;
    rms_residual  = sqrt(mean(residual.^2));

    if verbose
        fprintf('\n--- DEBUG: CFO ESTIMATION ---\n');
        fprintf('slope (rad/symbol) = %.6f\n', slope);
        if ~isempty(CFO_true)
            fprintf('CFO uoc luong      = %.4f Hz (CFO that = %.4f Hz)\n', cfo_est, CFO_true);
        else
            fprintf('CFO uoc luong      = %.4f Hz\n', cfo_est);
        end
        fprintf('RMS residual sau LS fit = %.6f rad\n', rms_residual);
    end

    N_payload      = length(rx_payload);
    n_idx          = (0:N_payload-1).';
    cfo_correction = exp(-1j * 2*pi * cfo_est * Tsym * n_idx);
    rx_payload_cfo = rx_payload .* cfo_correction;

    if do_plot
        figure;
        scatter(real(rx_payload_cfo), imag(rx_payload_cfo), 'b.'); hold on;
        scatter([1 -1 1 -1], [1 1 -1 -1], 'rx', 'LineWidth', 2);
        title('Chom sao sau khi bu CFO');
        xlabel('In-phase'); ylabel('Quadrature');
        grid on; axis square;
    end

    info.slope                  = slope;
    info.residual                = residual;
    info.rms_residual            = rms_residual;
    info.pilot_phase_raw         = pilot_phase_raw;
    info.pilot_phase_unwrapped   = pilot_phase_unwrapped;
    info.coeffs                  = coeffs;
end

function [rx_payload_eq, rx_data_eq, W_final, W_track, info] = rls_equalize_zc(rx_payload_cfo, pilot_indices, zc_pilots, data_indices, varargin)
    p = inputParser;
    addParameter(p, 'Lambda_warmup', 0.2);
    addParameter(p, 'Lambda_main', 0.85);
    addParameter(p, 'P_init', 100);
    addParameter(p, 'W_init', 1);
    addParameter(p, 'N_warmup', 50);
    addParameter(p, 'PlotConst', false);
    addParameter(p, 'Verbose', false);
    parse(p, varargin{:});

    lambda_wu   = p.Results.Lambda_warmup;
    lambda_main = p.Results.Lambda_main;
    P           = p.Results.P_init;
    W           = p.Results.W_init;
    n_warmup    = p.Results.N_warmup;
    do_plot     = p.Results.PlotConst;
    verbose     = p.Results.Verbose;

    rx_payload_cfo = rx_payload_cfo(:);
    pilot_indices  = pilot_indices(:);
    zc_pilots      = zc_pilots(:);
    data_indices   = data_indices(:);

    N = length(rx_payload_cfo);
    rx_payload_eq = zeros(N, 1);
    W_track       = zeros(N, 1);

    first_pilot_rx = rx_payload_cfo(pilot_indices(1));
    first_pilot_tx = zc_pilots(1);

    for it = 1:n_warmup
        y_wu = first_pilot_rx;
        d_wu = first_pilot_tx;

        y_mag_sq = real(y_wu)^2 + imag(y_wu)^2;
        k = (P * conj(y_wu)) / (lambda_wu + P * y_mag_sq);
        e = d_wu - W * y_wu;
        W = W + k * e;
        P = (P - k * y_wu * P) / lambda_wu;
    end

    if verbose
        fprintf('\n--- DEBUG: SAU WARM-UP ---\n');
        fprintf('W = %.4f%+.4fi, P = %.4f\n', real(W), imag(W), P);
    end

    pilot_cnt = 0;

    for i = 1:N
        y_k = rx_payload_cfo(i);

        if ismember(i, pilot_indices)
            pilot_cnt = pilot_cnt + 1;
            d_k = zc_pilots(pilot_cnt);
        else
            eq_temp = W * y_k;
            d_k = sign(real(eq_temp)) + 1j * sign(imag(eq_temp));
        end

        rx_payload_eq(i) = W * y_k;

        y_mag_sq = real(y_k)^2 + imag(y_k)^2;
        k = (P * conj(y_k)) / (lambda_main + P * y_mag_sq);
        e = d_k - rx_payload_eq(i);
        W = W + k * e;
        P = (P - k * y_k * P) / lambda_main;

        W_track(i) = W;
    end

    rx_data_eq = rx_payload_eq(data_indices);
    W_final    = W;

    if do_plot
        figure;
        scatter(real(rx_data_eq), imag(rx_data_eq), 'b.'); hold on;
        scatter([1 -1 1 -1], [1 1 -1 -1], 'rx', 'LineWidth', 2);
        title('Chom sao truoc xoay pha (Da bu CFO + can bang RLS bang ZC Pilot)');
        xlabel('In-phase'); ylabel('Quadrature');
        grid on; axis square;
    end

    info.P_final    = P;
    info.pilot_cnt  = pilot_cnt;
    info.lambda_wu  = lambda_wu;
    info.lambda_main = lambda_main;
end
