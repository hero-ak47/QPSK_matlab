%% THÔNG SỐ HỆ THỐNG
Fs = 48000;      % Tần số lấy mẫu (Hz)
Fc = 10000;      % Tần số sóng mang (Hz)
Rb = 400;        % Tốc độ bit (bps)
Rs = Rb/2;       % Tốc độ symbol cho QPSK (200 baud)
L  = Fs/Rs;      % Số mẫu trên một symbol (240 samples/symbol)

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

% 3. TẠO DỮ LIỆU NGẪU NHIÊN (400 symbols = 800 bits)
rng(123);                              % Seed cố định để RX tái tạo tính BER
num_data_sym = 400;
data_bits = randi([0 1], num_data_sym * 2, 1);
data_sym = (1 - 2*data_bits(1:2:end)) + 1j*(1 - 2*data_bits(2:2:end));

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
tx_symbols = [preamble_sym; guard; payload_sym;pss];

% 5. TẠO DÁNG XUNG VÀ ĐIỀU CHẾ LÊN BĂNG TẦN CƠ SỞ (BASEBAND)
beta = 0.5; % Hệ số Roll-off
span = 6;   % Độ dài bộ lọc
h_rrc = rcosdesign(beta, span, L, 'sqrt');

% Upsample và Lọc (Pulse Shaping)
tx_baseband = upfirdn(tx_symbols, h_rrc, L);

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
%sound(tx_signal_final, 48000);
disp('Đã phát xong!');

audiowrite('tx_signal_1.wav', tx_signal_final, 48000);
