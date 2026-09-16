%% THÔNG SỐ HỆ THỐNG
Fs = 48000;      % Tần số lấy mẫu (Hz)
Fc = 12000;      % Tần số sóng mang (Hz)
Rb = 400;        % Tốc độ bit (bps)
Rs = Rb/2;       % Tốc độ symbol cho QPSK (200 baud)
L  = Fs/Rs;      % Số mẫu trên một symbol (240 samples/symbol)

% ================= THÔNG SỐ MÃ HÓA KÊNH (MỚI THÊM) =================
% Convolutional code K=7, rate 1/2, generator chuẩn cong nghiep (171,133 octal)
constraint_length = 7;
code_generators   = [171 133];         % octal
trellis           = poly2trellis(constraint_length, code_generators);
n_tail            = constraint_length - 1;   % = 6 bit flush can them de dua trellis ve state 0

% So bit "cho chua" hien co: 400 symbol data * 2 bit/symbol = 800 bit
N_slots_available = 400 * 2;           % = 800

% Voi rate 1/2: coded_length = 2*(n_info + n_tail)
% => n_info = N_slots_available/2 - n_tail
n_info = N_slots_available/2 - n_tail;  % = 400 - 6 = 394 bit thong tin goc

fprintf('=== THONG SO MA HOA KENH ===\n');
fprintf('Constraint length K = %d, rate = 1/2, generator = [%d %d] (octal)\n', ...
        constraint_length, code_generators(1), code_generators(2));
fprintf('So bit thong tin goc (truoc encode): %d bit\n', n_info);
fprintf('So bit sau encode (bao gom %d bit flush): %d bit\n', n_tail, 2*(n_info+n_tail));
fprintf('Overhead bang thong: +100%% (rate 1/2)\n\n');

% 1. TẠO PREAMBLE ZADOFF-CHU (Đồng bộ khung)
Nzc = 63;
u   = 25;
n   = (0:Nzc-1).';
pss = exp(-1j*pi*u*n.*(n+1)/Nzc);
preamble_sym = pss;

% 2. TẠO CHUỖI PILOT ZADOFF-CHU (Đồng bộ kênh cho RLS)
N_pilots = 100;                        % Total 100 pilot symbols
u_p      = 7;                          % Root index cho Pilot ZC (phải trùng khớp với RX)
n_p      = (0:N_pilots-1).';
% Nhân với sqrt(2) để đồng bộ công suất trung bình với QPSK (±1 ±1j)
zc_pilots = sqrt(2) * exp(-1j * pi * u_p * n_p .* (n_p + 1) / N_pilots);

% ================= TẠO & MÃ HÓA DỮ LIỆU (ĐÃ SỬA) =================
rng(123);                              % Seed cố định để RX tái tạo tính BER

% ---- Tạo bit thông tin GỐC (trước khi encode, ngắn hơn trước) ----
% ============================================================
% TEXT -> BIT
% ============================================================

tx_text = 'RaDiMangNangLoiThe, Chua Thang Giac My Chua Ve BK';

% Chuyển từng ký tự ASCII thành 8 bit
tx_bytes = uint8(tx_text);

info_bits_full = reshape( ...
    de2bi(tx_bytes, 8, 'left-msb').', ...
    [], 1);

fprintf('TX text: %s\n', tx_text);
fprintf('So bit text: %d\n', length(info_bits_full));

% Pad 0 nếu thiếu, lỗi nếu vượt quá n_info
if length(info_bits_full) > n_info
    error('Text dai qua! Chi duoc toi da %d bit.', n_info);
end

info_bits = zeros(n_info, 1);
info_bits(1:length(info_bits_full)) = info_bits_full;

% Luu bit goc
fid = fopen('info_bits.bin', 'wb');
fwrite(fid, uint8(info_bits), 'uint8');
fclose(fid);

fprintf('Text -> %d bit\n', length(info_bits_full));
fprintf('Da pad thanh %d bit thong tin\n', n_info);

% ---- Thêm tail bits (flush) để đưa trellis về state 0 ----
info_bits_tail = [info_bits; zeros(n_tail, 1)];

% ---- Convolutional encode ----
coded_bits = convenc(info_bits_tail, trellis);   % length = 800, khop dung 400 symbol
fid = fopen('coded_bits.bin', 'wb');
fwrite(fid, uint8(coded_bits), 'uint8');
fclose(fid);

fprintf('Kiem tra: length(coded_bits) = %d (can = %d)\n', ...
        length(coded_bits), N_slots_available);
assert(length(coded_bits) == N_slots_available, ...
       'Loi: so bit sau encode khong khop voi so cho chua!');

% ================= INTERLEAVER (MỚI THÊM) =================
% Xao tron bit theo hoan vi ngau nhien CO DINH (biet truoc o ca TX/RX)
% -> pha vo tinh "cum" (burst) cua loi do Rayleigh fading
rng(2024);   % seed RIENG cho interleaver, khac seed cua du lieu (123)
interleaver_pattern = randperm(length(coded_bits));
coded_bits_il = coded_bits(interleaver_pattern);
fid = fopen('coded_bits_il.bin', 'wb');
fwrite(fid, uint8(coded_bits_il), 'uint8');
fclose(fid);

fprintf('Da ap dung interleaver (hoan vi ngau nhien %d bit)\n\n', length(coded_bits));

% ---- Map bit đã mã hóa VÀ ĐÃ INTERLEAVE sang QPSK ----
num_data_sym = 400;
data_sym = (1 - 2*coded_bits_il(1:2:end)) + 1j*(1 - 2*coded_bits_il(2:2:end));

% 4. LẮP RÁP PAYLOAD (1 Pilot ZC kèm 4 Data -> Khối 5 symbols, lặp lại 100 lần)
payload_sym = zeros(500, 1);
data_idx = 1;
for k = 1:100
    idx = (k-1)*5 + 1;
    payload_sym(idx)            = zc_pilots(k);                   % Chèn 1 symbol ZC pilot
    payload_sym(idx+1 : idx+4)  = data_sym(data_idx : data_idx+3);% Chèn 4 Data symbols
    data_idx = data_idx + 4;
end

% Tổng hợp toàn bộ symbols của khung (Preamble + Guard + Payload)
u_g = 41;                              % root khác preamble (u=25) và pilot (u_p=7)
n_g = (0:13).';
guard = exp(-1j*pi*u_g*n_g.*(n_g+1)/14);   % |guard| = 1, biên độ liên tục
tx_symbols = [preamble_sym; guard; payload_sym; pss];

% 5. TẠO DÁNG XUNG VÀ ĐIỀU CHẾ BASEBAND
beta = 0.5;
span = 6;

h_rrc = rcosdesign(beta, span, L, 'sqrt');

% ============================================================
% TX RRC
% ============================================================
tx_baseband = upfirdn(tx_symbols, h_rrc, L);

% ============================================================
% KÊNH AWGN
% ============================================================
SNR_dB = 1;     % Có thể đổi: 5, 8, 10, 12, 15, 20...

rng(999);        % seed cố định để tái lập kết quả

sig_power = mean(abs(tx_baseband).^2);
noise_power = sig_power / (10^(SNR_dB/10));

noise = sqrt(noise_power/2) * ...
        (randn(size(tx_baseband)) + ...
         1j*randn(size(tx_baseband)));

rx_channel = tx_baseband + noise;

fprintf('\n=== AWGN CHANNEL ===\n');
fprintf('SNR = %.2f dB\n', SNR_dB);
fprintf('Signal power = %.6e\n', sig_power);
fprintf('Noise power  = %.6e\n', noise_power);

% ============================================================
% RX MATCHED FILTER
% ============================================================
rx_mf = conv(rx_channel, h_rrc, 'same');

fprintf('Da qua RX matched filter\n');

%% TAO capture.bin SAU CA TX RRC + RX MATCHED FILTER

CAPTURE_SEC = 5;
N_CAPTURE   = Fs * CAPTURE_SEC;

x = rx_mf(:);

% Normalize
x = x / max(abs(x));

% Pad 5 giay
capture = complex(zeros(N_CAPTURE,1));

N_copy = min(length(x), N_CAPTURE);
capture(1:N_copy) = x(1:N_copy);

% Q1.31
scale = 2147483647;

I = int32(round(real(capture) * scale));
Q = int32(round(imag(capture) * scale));

% I/Q interleaved
iq = zeros(2*N_CAPTURE,1,'int32');
iq(1:2:end) = I;
iq(2:2:end) = Q;

fid = fopen('capture.bin','wb','ieee-le');
assert(fid ~= -1);

fwrite(fid, iq, 'int32');
fclose(fid);

fprintf('capture.bin da tao: TX RRC -> RX matched filter\n');
% 6. ĐIỀU CHẾ SÓNG MANG (PASSBAND)
t = (0:length(tx_baseband)-1)' / Fs;
tx_passband = real(tx_baseband .* exp(1j * 2 * pi * Fc * t));

% Chuẩn hóa biên độ tín hiệu (0.8 để tránh vỡ tiếng trên loa)
tx_passband = tx_passband / max(abs(tx_passband)) * 0.8;

% Thêm khoảng lặng ở đầu và cuối (0.5 giây)
silence = zeros(Fs * 0.5, 1);
tx_signal_final = [silence; tx_passband; silence];

% 7. PHÁT TÍN HIỆU RA LOA
disp('Đang phát tín hiệu ra loa...');
sound(tx_signal_final, 48000);
disp('Đã phát xong!');
audiowrite('tx_signal_2.wav', tx_signal_final, 48000);

% ================= LƯU LẠI THÔNG TIN CẦN CHO RX GIẢI MÃ =================
% RX se can: trellis, n_info, n_tail, va info_bits (de tinh BER dung tren
% BIT GOC thay vi bit da encode)
save('tx_coding_info.mat', 'trellis', 'n_info', 'n_tail', 'info_bits', ...
     'constraint_length', 'code_generators', 'interleaver_pattern');
fprintf('\nDa luu thong tin ma hoa vao tx_coding_info.mat (RX se can file nay de decode).\n');

