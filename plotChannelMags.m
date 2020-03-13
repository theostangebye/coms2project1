% Part II asks us to send data through bandlimited channels whose bw is the
% 10 dB bandwidth of the 10 channel datasets.

% Here, we plot Bode plots for the 10 channels to examine their 10 dB
% Bandwidths.

clear
close all                   % Close all current figure
load("comms432proj1.mat");  % Load channel data

for i=1:10
    subplot(2,5,i);
    plot(f,10*log10(abs(Cf(:,i))))
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'linear')
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.YAxis.TickLabelFormat = '%.2f';
    xlabel('f [MHz]');
    ylabel('|C(f)| [dB]');
    title(['Attenuation Length ',num2str(cz(i),4)]);
end