% Theo and Stephen

close all                   % Close all current figure
load("comms432proj1.mat");  % Load channel data
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 600000;                 % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor
rng default                 % Use default random number generator

% Constants defined by theo - should be changed?
P = 10e-3; % power is 10mWatts
No = 10e-14; % noise floor

sym_rates = f*1e6;              % Datarates from 0 to 1.1 Gb/s.  f given by clarkson dataset
bers = sym_rates;                % allowcate space for BER results (by copying sym_rates).
snrs = sym_rates;

for i=1:size(sym_rates,1)        % Iterate over the 2nd dimension of sym_rates.
    
    Ts = 1./sym_rates(i);         % symbol time
    snr_bit_db = 10*log10((P*Ts)/(k*No));   % snr per bit
    snr = snr_bit_db + 10*log10(k) - 10*log10(numSamplesPerSymbol); % converts bitwise ebno to symbolwise ebno
    
    dataIn = randi([0 1],n,1);  % Generate vector of random binary data
    
    % Fulfils case 1 of project instructions
    dataInMatrix = reshape(dataIn,length(dataIn)/k,k);  % Reshape data into binary k-tuples, k = log2(M)
    dataSymbolsIn = bi2de(dataInMatrix);                % Convert to integers
    dataMod = qammod(dataSymbolsIn,M,'bin');            % Binary coding, phase offset = 0
    receivedSignal = awgn(dataMod,snr,'measured');      % send data through awgn channel
    dataSymbolsOut = qamdemod(receivedSignal,M,'bin');  % detect symbols
    dataOutMatrix = de2bi(dataSymbolsOut,k);            % symbols to binary
    dataOut = dataOutMatrix(:);                         % Return data in column vector
    [numErrors,ber] = biterr(dataIn,dataOut);           % compute BER
    bers(i) = ber;                                      % save BER value
    snrs(i) = snr;
    
end

% make a nice plot below:
plot(snrs,bers)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'linear')
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
xlabel('SNR');
ylabel('BER');
hold on
title('BER vs. SNR');
% HOLD ON, CHANGE M, run again and then do hold off.
legend('16 QAM')