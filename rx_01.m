
clc; clear; close all
% SECTION 1: BẮT ĐẦU THU TÍN HIỆU
% disp('Đang ghi âm... Hãy phát tín hiệu ở máy kia.');
% r = audiorecorder(48000, 24, 1);
% record(r);

 %% SECTION 2: DỪNG THU VÀ XỬ LÝ TÍN HIỆU
% stop(r);
% disp('Đã dừng ghi âm. Đang xử lý tín hiệu...');
%y = getaudiodata(r, 'double');
% 1. THÔNG SỐ VÀ TÁI TẠO PREAMBLE LÀM CHUẨN MẪU
Fs = 48000; Fc = 10000; 
Rs = 200;   % bit rate
L = Fs/Rs;  % numbers of sample in 1s

% RRC
beta = 0.5; span = 6;
h_rrc = rcosdesign(beta, span, L, 'sqrt');

% Zadoff-Chu
Nzc = 63;
u   = 25;
Lsym = Nzc * L;              % Số mẫu của MỘT khối preamble ZC sau pulse shaping
n   = (0:Nzc-1).';
pss = exp(-1j*pi*u*n.*(n+1)/Nzc);

% preable laf ZC de dong bo
preamble_sym_rx = pss;
preamble_bb_ideal = upfirdn(preamble_sym_rx, h_rrc, L);

% ================= NẠP TÍN HIỆU =================
%load('tx_signal_2.mat', 'y');
%-------------------------------------
[y, Fs_read] = audioread('tx_signal_1.wav');
if Fs_read ~= Fs
    warning('Fs trong file (%d) khác Fs khai báo (%d)!', Fs_read, Fs);
end
disp('Đã nạp tín hiệu từ tx_signal_1.wav');

% Khu offset
y = y - mean(y);

% a_true = 0.005;   % giả lập Doppler wideband thật (a = v/c)
% t_orig = (0:length(y)-1)'/Fs;
% t_warp = t_orig*(1+a_true);
% y = interp1(t_orig, y, t_warp, 'spline', 2);

rayleighChan = comm.RayleighChannel( ...
    'SampleRate', 48000, ...
    'PathDelays', [0 0.001 0.002], ...
    'AveragePathGains', [0 -9 -10], ...
    'MaximumDopplerShift', 1.5);
y = rayleighChan(y);

% mo phong kenh truyen : CFO 
CFO_true = 0.5; 
y = y .* exp(1j*2*pi*CFO_true*(0:length(y)-1)'/Fs);
%awgn
SNR_dB = 25;
y = awgn(y, SNR_dB);
plot(y);



% 2. TIỀN XỬ LÝ: LỌC BĂNG THÔNG & HẠ TẦN XUỐNG BASEBAND
[b, a] = butter(5, [(Fc - 1000) (Fc + 1000)]/(Fs/2), 'bandpass');
y_filt = filter(b, a, y);

t_rx = (0:length(y_filt)-1)' / Fs;
y_bb = y_filt .* exp(-1j * 2 * pi * Fc * t_rx);

y_mf = conv(y_bb, h_rrc, 'same');

% 3. ĐỒNG BỘ KHUNG (FRAME SYNCHRONIZATION) — PHIÊN BẢN AN TOÀN
Nzc       = 63;
N_guard   = 14;
N_payload = 500;
N_between_symbols = Nzc + N_guard + N_payload;   % = 577 symbols
T0 = N_between_symbols * L / Fs;
T0_samples = N_between_symbols * L;

[xc, lags] = xcorr(y_mf, preamble_bb_ideal);
xc_abs = abs(xc);
noise_floor = median(xc_abs);

% ---- Tìm ĐỈNH 1 (preamble đầu) — CHỈ trong nửa đầu tín hiệu ----
search_limit_lag = round(length(y_mf) / 2);
valid_idx1 = find(lags >= 0 & lags <= search_limit_lag);
if isempty(valid_idx1)
    error('Vùng tìm kiếm đỉnh 1 không hợp lệ -> kiểm tra lại độ dài tín hiệu!');
end

[peak1_val, rel_idx1] = max(xc_abs(valid_idx1));
idx1 = valid_idx1(rel_idx1);
lag1 = lags(idx1);
start_sample = lag1 + 1;

fprintf('\n--- DEBUG: FRAME SYNC (đỉnh 1, nửa đầu tín hiệu) ---\n');
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
%------------------------------------------------------------------
fprintf('\n--- DEBUG: 2-PILOT DOPPLER ESTIMATION ---\n');
fprintf('T0 (khoảng cách danh định giữa 2 preamble) = %.6f s\n', T0);

% Vùng tìm đỉnh 2: bắt đầu từ NGAY SAU vùng đã dùng cho đỉnh 1,
% và không vượt quá độ dài tín hiệu
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

% ---- Tính hệ số Doppler ----
sample_diff = lag2 - lag1;
T_rx = sample_diff / Fs;
a_hat_2pilot = T0/T_rx - 1;

fprintf('T_rx (đo được)   = %.6f s\n', T_rx);
fprintf('a_hat (2-pilot)  = %.6f\n', a_hat_2pilot);
fprintf('f_D tương đương ~ %.3f Hz (tại Fc = %d Hz)\n', a_hat_2pilot*Fc, Fc);
fprintf('Vận tốc tương đương ~ %.3f m/s (c = 1500 m/s)\n', a_hat_2pilot*1500);
%------------------------------------------------------------------
% ---- ÁP DỤNG BÙ DOPPLER TIME-SCALE ----
y = doppler_resample(y, a_hat_2pilot, Fs);
fprintf('\nĐã bù Doppler time-scale, độ dài tín hiệu sau bù: %d mẫu\n', length(y));

% ---- LẶP LẠI TOÀN BỘ TIỀN XỬ LÝ + ĐỒNG BỘ KHUNG TRÊN TÍN HIỆU ĐÃ BÙ ----
y_filt = filter(b, a, y);
t_rx = (0:length(y_filt)-1)' / Fs;
y_bb = y_filt .* exp(-1j * 2 * pi * Fc * t_rx);
y_mf = conv(y_bb, h_rrc, 'same');

[xc, lags] = xcorr(y_mf, preamble_bb_ideal);
xc_abs = abs(xc);
noise_floor = median(xc_abs);

search_limit_lag = round(length(y_mf) / 2);
valid_idx1 = find(lags >= 0 & lags <= search_limit_lag);
[peak1_val, rel_idx1] = max(xc_abs(valid_idx1));
idx1 = valid_idx1(rel_idx1);
start_sample = lags(idx1) + 1;   % <-- start_sample MỚI, trên tín hiệu ĐÃ BÙ

fprintf('\n--- FRAME SYNC LẠI sau khi bù Doppler ---\n');
fprintf('start_sample (mới) = %d\n', start_sample);
fprintf('Peak1/Noise = %.2f\n', peak1_val/noise_floor);

if start_sample < 1 || (start_sample + (Nzc+N_guard+N_payload)*L) > length(y_mf)
    error('Không tìm thấy khung tín hiệu nguyên vẹn sau khi bù Doppler!');
end


% 4. TRÍCH XUẤT SYMBOL (DOWNSAMPLING) -- lấy mẫu đại diện cho các symbol IQ
total_sym_len = 577 + 63;
skip = 77;
rx_symbols = zeros(total_sym_len, 1);

first_sym_offset = (span*L/2) + 1;  % lấy mẫu ở chính giữa symbol

for k = 1:total_sym_len
    idx = start_sample + first_sym_offset + (k-1)*L - 1;
    rx_symbols(k) = y_mf(idx);
end

rx_payload = rx_symbols(skip+1:end - Nzc); % Bỏ qua preamble + guard  -> thu được các symbol pilot + data

% ---- Vị trí pilot & chuỗi pilot ZC (khai báo TRƯỚC để dùng cho cả CFO và RLS) ----
pilot_indices = 1:5:500;             % 100 vị trí pilot
data_indices  = setdiff(1:500, pilot_indices); % hàm setdiff trả về hiệu của 2 tập hợp -> index của data là các số còn lại từ 1 -> 500

N_pilots = length(pilot_indices);    % = 100 pilot symbols
u_p = 7;
n_p = (0:N_pilots-1).';
zc_pilots = sqrt(2) * exp(-1j * pi * u_p * n_p .* (n_p + 1) / N_pilots);

%-------------------------------------------------------------------------------------------------------------------------------
% 5. ƯỚC LƯỢNG & BÙ CFO BẰNG LEAST-SQUARES TRÊN PHA PILOT ZC
%-------------------------------------------------------------------------------------------------------------------------------
[cfo_est, rx_payload_cfo, info] = estimate_and_correct_cfo(...
    rx_payload, pilot_indices, zc_pilots, L, Fs, ...
    'CFO_true', CFO_true, 'Verbose', true, 'PlotConst', true);
%rx_payload_cfo = rx_payload;

%----------------------------------------------------------------
% 6. CÂN BẰNG KÊNH BẰNG RLS (DÙNG PILOT ZC) - CHẠY TRÊN rx_payload_cfo
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

% ---- Kiểm tra hội tụ: vẽ |W| và pha theo symbol index ----
figure;
subplot(2,1,1);
plot(abs(W_track), 'LineWidth', 1.5); grid on;
title('Bien do trong so |W| cap nhat qua 500 symbol (sau bu CFO)');
xlabel('Symbol index'); ylabel('|W|');

subplot(2,1,2);
plot(rad2deg(angle(W_track)), 'LineWidth', 1.5); grid on;
title('Pha cua W cap nhat qua 500 symbol (sau bu CFO)');
xlabel('Symbol index'); ylabel('Pha (do)');

% 7. GIẢI ĐIỀU CHẾ QPSK (DEMAPPING) VÀ SO SÁNH
rx_bits = zeros(length(rx_data_eq)*2, 1);
rx_bits(1:2:end) = real(rx_data_eq) < 0;
rx_bits(2:2:end) = imag(rx_data_eq) < 0;

rng(123);
tx_data_bits_ideal = randi([0 1], 400 * 2, 1);

num_errors = sum(rx_bits ~= tx_data_bits_ideal);
BER = num_errors / length(tx_data_bits_ideal);

% 8. HIỂN THỊ KẾT QUẢ
disp('--- KẾT QUẢ TRUYỀN NHẬN (Bu CFO (LS pilot ZC) + Can bang RLS) ---');
fprintf('Số bits truyền: %d\n', length(tx_data_bits_ideal));
fprintf('Số bits lỗi: %d\n', num_errors);
fprintf('Tỉ lệ lỗi bit (BER): %f\n', BER);

figure;
scatter(real(rx_data_eq), imag(rx_data_eq), 'b.'); hold on;
scatter([1 -1 1 -1], [1 1 -1 -1], 'rx', 'LineWidth', 2);
title('Chòm sao tín hiệu nhận được (Bù CFO + Cân bằng RLS bằng ZC Pilot)');
xlabel('In-phase'); ylabel('Quadrature');
grid on; axis square;





% Tạo struct chứa metadata
% info = struct();
% info.description   = 'Tin hieu QPSK preamble ZC, Fs thiet ke = 48000, Fs header ghi sai 48050 (co tinh)';
% info.Fs_design      = 48000;
% info.Fs_header      = 48100;
% info.date_created   = datetime('now');
% info.Fc             = 8000;          % Doppler thật (nếu có)
% info.SNR_dB         = 10;
% info.notes          = 'Dung de test doppler_resample, KHONG qua loa-mic that';
% 
% % Lưu cả y lẫn info vào cùng 1 file .mat
% save('tx_signal_6.mat', 'y', 'info');

function y_out = doppler_resample(y, a_try, Fs, fc)
    % y    : Tín hiệu đầu vào (Real hoặc Complex IQ)
    % a_try: Hệ số Doppler (a = v/c)
    % Fs   : Tần số lấy mẫu (Hz)
    % fc   : Tần số sóng mang (Hz). Bỏ trống hoặc đặt = 0 nếu y là Passband.

    if nargin < 4 || isempty(fc)
        fc = 0; % Mặc định xem như tín hiệu passband hoặc đã hạ băng hoàn hảo
    end

    N = length(y);
    t_orig = (0:N-1)' / Fs;        % Trục thời gian gốc

    % 1. Đổi trục thời gian chính xác: t_new = t / (1 + a)
    t_new = t_orig / (1 + a_try);

    % 2. Nội suy và cho phép ngoại suy ('extrap') tránh bị tràn số 0 ở biên
    y_resampled = interp1(t_orig, y, t_new, 'spline', 'extrap');

    % 3. Bù xoay pha sóng mang (Carrier Phase Shift) nếu là tín hiệu baseband
    if fc ~= 0
        y_out = y_resampled .* exp(-1i * 2 * pi * fc * a_try * t_orig);
    else
        y_out = y_resampled;
    end
end

function [cfo_est, rx_payload_cfo, info] = estimate_and_correct_cfo(rx_payload, pilot_indices, zc_pilots, L, Fs, varargin)
% ESTIMATE_AND_CORRECT_CFO  Uoc luong CFO tu pha cac pilot Zadoff-Chu (ZC)
% bang LS fit tuyen tinh, sau do bu CFO cho toan bo payload.
%
% INPUT:
%   rx_payload    - vector cot (Nx1) cac symbol payload da thu duoc (IQ)
%   pilot_indices - vi tri (1-based) cua cac pilot trong rx_payload
%   zc_pilots     - chuoi ZC ly tuong tuong ung voi pilot_indices
%   L             - do dai (so mau) cua 1 symbol -> dung tinh Tsym = L/Fs
%   Fs            - tan so lay mau (Hz)
%
% OPTIONAL NAME-VALUE:
%   'CFO_true'  - gia tri CFO that (Hz), chi de in debug/so sanh (default: [])
%   'Verbose'   - true/false, co in debug ra console khong (default: true)
%   'PlotConst' - true/false, co ve chom sao sau bu CFO khong (default: false)
%
% OUTPUT:
%   cfo_est        - CFO uoc luong (Hz)
%   rx_payload_cfo - rx_payload sau khi da bu CFO
%   info           - struct chua thong tin debug:
%                     .slope, .residual, .rms_residual,
%                     .pilot_phase_raw, .pilot_phase_unwrapped,
%                     .coeffs

    % ---- parse optional args ----
    p = inputParser;
    addParameter(p, 'CFO_true', []);
    addParameter(p, 'Verbose', true);
    addParameter(p, 'PlotConst', false);
    parse(p, varargin{:});
    CFO_true  = p.Results.CFO_true;
    verbose   = p.Results.Verbose;
    do_plot   = p.Results.PlotConst;

    % ---- dam bao vector cot ----
    rx_payload    = rx_payload(:);
    pilot_indices = pilot_indices(:);
    zc_pilots     = zc_pilots(:);

    % ---- lay pilot thu duoc tai cac vi tri pilot_indices ----
    rx_pilots_raw = rx_payload(pilot_indices);

    % ---- pha cua pilot thu duoc so voi pilot ZC ly tuong ----
    pilot_phase_raw = angle(rx_pilots_raw ./ zc_pilots);

    if verbose
        fprintf('\n--- DEBUG: PILOT PHASE (TRUOC UNWRAP) ---\n');
        disp(pilot_phase_raw(1:min(10, end)).');
    end

    % ---- unwrap pha de bu k*2*pi ----
    pilot_phase_unwrapped = unwrap(pilot_phase_raw);

    % ---- ma tran thiet ke A = [1 n], phi_CFO = A*[a;b] ----
    idx_col = pilot_indices - 1;               % 0-based index lam "thoi gian"
    A = [ones(length(idx_col), 1), idx_col];

    % ---- nghiem LS: phi = (A'A)^-1 A' phi_unwrap ----
    coeffs = A \ pilot_phase_unwrapped;
    slope  = coeffs(2);

    % ---- tinh CFO tu do doc (slope) ----
    Tsym    = L / Fs;
    cfo_est = slope / (2*pi*Tsym);

    % ---- residual sau fit ----
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

    % ---- bu CFO cho toan bo payload ----
    N_payload      = length(rx_payload);
    n_idx          = (0:N_payload-1).';
    cfo_correction = exp(-1j * 2*pi * cfo_est * Tsym * n_idx);
    rx_payload_cfo = rx_payload .* cfo_correction;

    % ---- ve chom sao (tuy chon) ----
    if do_plot
        figure;
        scatter(real(rx_payload_cfo), imag(rx_payload_cfo), 'b.'); hold on;
        scatter([1 -1 1 -1], [1 1 -1 -1], 'rx', 'LineWidth', 2);
        title('Chom sao sau khi bu CFO');
        xlabel('In-phase'); ylabel('Quadrature');
        grid on; axis square;
    end

    % ---- dong goi info de debug them neu can ----
    info.slope                  = slope;
    info.residual                = residual;
    info.rms_residual            = rms_residual;
    info.pilot_phase_raw         = pilot_phase_raw;
    info.pilot_phase_unwrapped   = pilot_phase_unwrapped;
    info.coeffs                  = coeffs;

end

function [rx_payload_eq, rx_data_eq, W_final, W_track, info] = rls_equalize_zc(rx_payload_cfo, pilot_indices, zc_pilots, data_indices, varargin)
% RLS_EQUALIZE_ZC  Can bang kenh 1-tap bang RLS (Recursive Least Squares),
% dung Decision-Directed cho symbol data va ZC pilot lam tham chieu that
% cho symbol pilot. Co giai doan warm-up tren pilot dau tien de hoi tu nhanh.
%
% INPUT:
%   rx_payload_cfo - vector cot (Nx1) tin hieu payload DA BU CFO
%   pilot_indices  - vi tri (1-based) cac pilot trong rx_payload_cfo
%   zc_pilots      - chuoi ZC ly tuong tuong ung voi pilot_indices
%   data_indices   - vi tri (1-based) cac symbol data trong rx_payload_cfo
%
% OPTIONAL NAME-VALUE:
%   'Lambda_warmup' - he so quen dung trong vong warm-up (default: 0.2)
%   'Lambda_main'   - he so quen dung trong vong RLS chinh (default: 0.65)
%   'P_init'        - gia tri khoi tao P (default: 100)
%   'W_init'        - gia tri khoi tao trong so W (default: 1)
%   'N_warmup'      - so lan lap warm-up tren pilot dau tien (default: 50)
%   'PlotConst'      - true/false, ve chom sao rx_data_eq sau cung (default: false)
%   'Verbose'        - true/false, in thong tin debug (default: false)
%
% OUTPUT:
%   rx_payload_eq - toan bo payload sau khi can bang RLS (Nx1)
%   rx_data_eq    - cac symbol data (tai data_indices) sau can bang
%   W_final       - trong so RLS cuoi cung sau khi chay het payload
%   W_track       - lich su trong so W qua tung symbol (Nx1)
%   info          - struct chua P cuoi cung, so pilot da dung, v.v.

    % ---- parse optional args ----
    p = inputParser;
    addParameter(p, 'Lambda_warmup', 0.2);
    addParameter(p, 'Lambda_main', 0.65);
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

    % ---- dam bao vector cot ----
    rx_payload_cfo = rx_payload_cfo(:);
    pilot_indices  = pilot_indices(:);
    zc_pilots      = zc_pilots(:);
    data_indices   = data_indices(:);

    N = length(rx_payload_cfo);
    rx_payload_eq = zeros(N, 1);
    W_track       = zeros(N, 1);

    % ==== Warm-up tai pilot ZC dau tien de hoi tu nhanh ====
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

    % ==== RLS lien tuc qua toan bo payload da bu CFO ====
    % (Decision-Directed cho data, ZC Pilot lam tham chieu that cho pilot)
    pilot_cnt = 0;

    for i = 1:N
        y_k = rx_payload_cfo(i);   % tin hieu DA BU CFO

        % 1. Lay du lieu tham chieu (Reference)
        if ismember(i, pilot_indices)      % neu dang xet symbol pilot
            pilot_cnt = pilot_cnt + 1;
            d_k = zc_pilots(pilot_cnt);    % pilot tham chieu da biet
        else                                % neu dang xet symbol data
            eq_temp = W * y_k;
            d_k = sign(real(eq_temp)) + 1j * sign(imag(eq_temp)); % symbol tham chieu du doan
        end

        % 2. Can bang tin hieu nhan duoc, luu vao mang ket qua
        rx_payload_eq(i) = W * y_k;

        % 3. Cap nhat RLS cho symbol tiep theo
        y_mag_sq = real(y_k)^2 + imag(y_k)^2;

        % do loi Kalman: quyet dinh tin sai so e den muc nao
        k = (P * conj(y_k)) / (lambda_main + P * y_mag_sq);

        % sai so uoc luong
        e = d_k - rx_payload_eq(i);

        % cap nhat trong so RLS moi
        W = W + k * e;

        % cap nhat P: P cang lon, thuat toan cang tu tin cap nhat manh
        % -> hoi tu nhanh nhung de dao dong
        P = (P - k * y_k * P) / lambda_main;

        W_track(i) = W;    % luu lai trong so de plot
    end

    rx_data_eq = rx_payload_eq(data_indices);
    W_final    = W;

    % ---- ve chom sao (tuy chon) ----
    if do_plot
        figure;
        scatter(real(rx_data_eq), imag(rx_data_eq), 'b.'); hold on;
        scatter([1 -1 1 -1], [1 1 -1 -1], 'rx', 'LineWidth', 2);
        title('Chom sao truoc xoay pha (Da bu CFO + can bang RLS bang ZC Pilot)');
        xlabel('In-phase'); ylabel('Quadrature');
        grid on; axis square;
    end

    % ---- dong goi info de debug them neu can ----
    info.P_final    = P;
    info.pilot_cnt  = pilot_cnt;
    info.lambda_wu  = lambda_wu;
    info.lambda_main = lambda_main;

end




