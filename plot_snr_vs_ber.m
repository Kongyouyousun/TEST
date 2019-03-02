%% snr 10:20
figure(1)
% semilogy(SNRIn(1:end),ber_snr(1:end), 'bo-','LineWidth',2);
semilogy(SNRIn(1:end),ber_snr(1:end), 'k');
grid on
xlabel(' SNR [dB]','FontSize',12)
  ylabel('bit error rate (BER)','FontSize',12)
axis([10 15 1e-6 1*1e-1])
% axis tight
