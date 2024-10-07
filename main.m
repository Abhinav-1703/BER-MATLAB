% Parameters
N = 10^6; % Number of bits to transmit
Eb_N0_dB = 0:10; % Range of Eb/N0 values in dB

% Generate random binary data
data = randi([0 1], N, 1); 

% Initialize BER storage for all modulation types
BER_BASK = zeros(length(Eb_N0_dB),1);
BER_BPSK = zeros(length(Eb_N0_dB),1);
BER_BFSK = zeros(length(Eb_N0_dB),1);
BER_QPSK = zeros(length(Eb_N0_dB),1);

for i = 1:length(Eb_N0_dB)
    % Convert Eb/N0 from dB to linear scale
    Eb_N0 = 10^(Eb_N0_dB(i)/10); 
    noise_std_dev = 1/sqrt(2*Eb_N0); 
    
    %% BASK Modulation and Demodulation
    tx_BASK = data; % BASK: 0 -> 0, 1 -> 1
    noise_BASK = noise_std_dev * randn(N, 1);
    rx_BASK = tx_BASK + noise_BASK; 
    demod_BASK = rx_BASK >= 0.5; % Decision threshold = 0.5
    BER_BASK(i) = sum(data ~= demod_BASK) / N;
    
    %% BPSK Modulation and Demodulation
    tx_BPSK = 2*data - 1; % BPSK: 0 -> -1, 1 -> +1
    noise_BPSK = noise_std_dev * randn(N, 1);
    rx_BPSK = tx_BPSK + noise_BPSK; 
    demod_BPSK = rx_BPSK >= 0; % Decision threshold = 0
    BER_BPSK(i) = sum(data ~= demod_BPSK) / N;
    
    %% BFSK Modulation and Demodulation
    tx_BFSK = 2*data - 1; % BFSK: 0 -> -1, 1 -> +1 (frequency representation simplified)
    noise_BFSK = noise_std_dev * randn(N, 1);
    rx_BFSK = tx_BFSK + noise_BFSK;
    demod_BFSK = rx_BFSK >= 0; % Decision threshold = 0
    BER_BFSK(i) = sum(data ~= demod_BFSK) / N;
    
    %% QPSK Modulation and Demodulation
    data_QPSK = reshape(data, [], 2); % Reshape to pairs of bits
    tx_QPSK = 1/sqrt(2) * (2*data_QPSK(:,1) - 1 + 1i * (2*data_QPSK(:,2) - 1)); % QPSK Modulation
    noise_QPSK = (noise_std_dev/sqrt(2)) * (randn(length(tx_QPSK),1) + 1i*randn(length(tx_QPSK),1));
    rx_QPSK = tx_QPSK + noise_QPSK;
    demod_QPSK_real = real(rx_QPSK) >= 0; % Decision for real part
    demod_QPSK_imag = imag(rx_QPSK) >= 0; % Decision for imaginary part
    demod_QPSK = [demod_QPSK_real demod_QPSK_imag];
    demod_QPSK = demod_QPSK(:); % Reshape to original bit vector
    BER_QPSK(i) = sum(data ~= demod_QPSK) / N;
end

% Create subplots for each modulation scheme
figure;

% BASK
subplot(2, 2, 1);
semilogy(Eb_N0_dB, BER_BASK, 'r-o');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs E_b/N_0 for BASK');
grid on;

% BPSK
subplot(2, 2, 2);
semilogy(Eb_N0_dB, BER_BPSK, 'b-s');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs E_b/N_0 for BPSK');
grid on;

% BFSK
subplot(2, 2, 3);
semilogy(Eb_N0_dB, BER_BFSK, 'g-^');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs E_b/N_0 for BFSK');
grid on;

% QPSK
subplot(2, 2, 4);
semilogy(Eb_N0_dB, BER_QPSK, 'k-d');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs E_b/N_0 for QPSK');
grid on;
