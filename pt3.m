% Theo and Stephen
clear
hold off
close all                   % Close all current figure
load("comms432proj1.mat");  % Load channel data
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 100000;                 % Number of bits to process
numSamplesPerSymbol = 1;    % Oversampling factor
span = 10;
rolloff = 0.25;
nSamp = 1;
rng default                 % Use default random number generator

styles = {'-.b*','--ko',':rs','-.g*','--ro',':ms','-.m*','--co',':bs','-.k*'};

% Constants defined by theo - should be changed?
P = 10e-4; % power is 10mWatts
No = 10e-14; % noise floor

% sym_rates = f*1e6;        % Datarates from 0 to 1.1 Gb/s.  f given by clarkson dataset
snrs = 1:30;
bers = snrs;                % allowcate space for BER results (by copying sym_rates).
losses = cz;


% TODO: REPLACE THIS CODE WITH CODE FROM HERE: https://www.mathworks.com/help/comm/gs/use-pulse-shaping-on-16-qam-signal.html
% Pulse Shaping Filters:
txfilter = comm.RaisedCosineTransmitFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'OutputSamplesPerSymbol',nSamp);
rxfilter = comm.RaisedCosineReceiveFilter('RolloffFactor',rolloff, ...
    'FilterSpanInSymbols',span,'InputSamplesPerSymbol',nSamp, ...
    'DecimationFactor',nSamp);

dataIn = randi([0 1],n,1);  % Generate vector of random binary data
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);  % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInMatrix);                % Convert to integers
dataMod = qammod(dataSymbolsIn,M,'bin');            % Binary coding, phase offset = 0

dataMod = txfilter(dataMod);

% Prep frequency data for fft.
f = f'.*1e6;                 % Convert from GHz to Hz.
dif = f(2)-f(1);            % distance between measurements.
f2 = [0:dif:f(1)-dif,f];    % create data for beginning of vector.

% 10 non-distortingband limited channels... (see handout)
for nn=1:size(cz,2) % iterate over each of 10 turbidity levels
        
    ch = Cf(:,nn)';    % Copy channel data for this sim run.
    ch = ch * (1/mean(abs(ch)));
     z_pad = zeros(1,size(f2,2) - size(ch,2));
     ch = [z_pad,ch]; 
     channel_t = ifft(ch, 'symmetric');


    for i=1:length(snrs)        % Iterate over the 2nd dimension of sym_rates.

        snr = snrs(i);
        
        receivedSignal = filter(channel_t,1,dataMod);
        
        receivedSignal = awgn(receivedSignal, snr,'measured');      % send data through awgn channel
        
        receivedSignal = rxfilter(receivedSignal);
        
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
    title('BER vs. SNR for Different Band-limited Distorting Channels.');
    legend(split(num2str(cz)))
    
    hold off;
%     figure
    
    fvtool(channel_t,'Analysis','freq');                                                  % Plot Mag of channel
    fvtool(channel_t,'Analysis','phase','Analysis','freq');                               % Plot phase  channel
    