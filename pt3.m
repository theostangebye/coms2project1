% Theo and Stephen
clear
hold off
close all                   % Close all current figure
load("comms432proj1.mat");  % Load channel data
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 100000;                 % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor
rng default                 % Use default random number generator

styles = {'-.b*','--ko',':rs','-.g*','--ro',':ms','-.m*','--co',':bs','-.k*'};

% Constants defined by theo - should be changed?
P = 10e-4; % power is 10mWatts
No = 10e-14; % noise floor

% sym_rates = f*1e6;        % Datarates from 0 to 1.1 Gb/s.  f given by clarkson dataset
snrs = 1:30;
bers = snrs;                % allowcate space for BER results (by copying sym_rates).
losses = cz;

dataIn = randi([0 1],n,1);  % Generate vector of random binary data
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);  % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInMatrix);                % Convert to integers
dataMod = qammod(dataSymbolsIn,M,'bin');            % Binary coding, phase offset = 0

% 10 non-distortingband limited channels... (see handout)
for nn=1:size(cz,2) % iterate over each of 10 turbidity levels
    
    % DISTORTING FILTER DESIGN
    dstrt_chnl = fdesign.arbmagnphase('N,F,H',100,f2./(max(f2)),Cf(:,nn)'); % estimate channel
    dstrt_fltr = design(dstrt_chnl);          % construct filter to emulate estimated channel
    fvtool(dstrt_fltr);                                                 % Plot Mag of channel
    fvtool(dstrt_fltr,'Analysis','phase');                              % Plot phase  channel
    % distort data like this: 
    % filter(dstrt_fltr, dataMod)
    
    for i=1:length(snrs)        % Iterate over the 2nd dimension of sym_rates.

        loss = Wc(nn)*1e6/1100e6;
        lossDB = 10*log10(loss);
        losses(nn) = lossDB;

        snr = snrs(i);
        
        receivedSignal = awgn(dataMod,snr + lossDB,'measured');      % send data through awgn channel
        dataSymbolsOut = qamdemod(receivedSignal,M,'bin');  % detect symbols
        dataOutMatrix = de2bi(dataSymbolsOut,k);            % symbols to binary
        dataOut = dataOutMatrix(:);                         % Return data in column vector
        [numErrors,ber] = biterr(dataIn,dataOut);           % compute BER
        bers(i) = ber;                                      % save BER value
        snrs(i) = snr; 
    end
    
    plot(snrs,bers,styles{nn})
    hold on
end

    % make a nice plot below:
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'linear')
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    xlabel('SNR');
    ylabel('BER');
    hold on
    title('BER vs. SNR for Different Band-limited Channel Losses.');
    legend(split(num2str(losses)))