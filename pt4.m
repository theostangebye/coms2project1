% Theo and Stephen
clear
hold off
close all                   % Close all current figure
load("comms432proj1.mat");  % Load channel data
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
n = 100000;                 % Number of bits to process
rng default                 % Use default random number generator

sps = 4; % Number of samples per symbol (oversampling factor)
filtlen = 10; % Filter length in symbols
rolloff = 0.25; % Filter rolloff factor

styles = {'-.b*','--ko',':rs','-.g*','--ro',':ms','-.m*','--co',':bs','-.k*'};

% Constants defined by theo - should be changed?
P = 10e-4; % power is 10mWatts
No = 10e-14; % noise floor

% sym_rates = f*1e6;        % Datarates from 0 to 1.1 Gb/s.  f given by clarkson dataset
snrs = 1:30;
bers = snrs;                % allowcate space for BER results (by copying sym_rates).
losses = cz;

Cf = Cf .* 1/mean(mean(abs(Cf)));

rrcFilter = rcosdesign(rolloff,filtlen,sps);

dataIn = randi([0 1],n,1);  % Generate vector of random binary data
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);  % Reshape data into binary k-tuples, k = log2(M)
dataSymbolsIn = bi2de(dataInMatrix);                % Convert to integers
dataMod = qammod(dataSymbolsIn,M,'bin');            % Binary coding, phase offset = 0

dataMod = upfirdn(dataMod,rrcFilter,sps,1);

% 10 non-distortingband limited channels... (see handout)
for nn=1:size(cz,2) % iterate over each of 10 turbidity levels
    
    % Determine loss of channel
    m = size(Cf,1); % number of data points per channel
    Chnl = Cf(:,nn); % channel data
    rrc_mag = freqz(rrcFilter,1.1e9,1024); % get mag and phase of rrc filter. m points.
    sum = 0;
    for j=1:m
        sum = sum + abs(rrc_mag(j))/(abs(Chnl(j)).^2);
    end
    sum = sum * 2; % account for negative frequencies.
    loss_rrc = 10*log10(sum);
    
    losses(nn) = loss_rrc;
    
    for i=1:length(snrs)        % Iterate over the 2nd dimension of sym_rates.

        snr = snrs(i);
         
        receivedSignal = awgn(dataMod, snr-loss_rrc,'measured');      % send data through awgn channel
        
        receivedSignal = upfirdn(receivedSignal,rrcFilter,1,sps); % Downsample and filter
        receivedSignal = receivedSignal(filtlen + 1:end - filtlen); % Account for delay
        
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
%     fvtool(dstrt_fltr);                                                  % Plot Mag of channel
%     fvtool(dstrt_fltr,'Analysis','phase');                               % Plot phase  channel
%     fvtool(eq_fltr);                                                  % Plot Mag of channel
%     fvtool(eq_fltr,'Analysis','phase');                               % Plot phase  channel
%     