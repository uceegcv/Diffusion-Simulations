%% Plot SNR graphs

xvalues = [20;30;40;50;60;100];


newdata = load('ave_error_5_3to200wavelengths_SNR50_origDPF_perfectE_SNRvals_concs2 .mat');
data_to_analyse = newdata.ave_error_SNR;
std_to_analyse = newdata.stand_dev_error_SNR;

    change_wavs_HbO_zero_error = data_to_analyse(:,1);
    change_wavs_HbR_zero_error = data_to_analyse(:,2);
    change_wavs_CCO_zero_error = data_to_analyse(:,3);
    std_wavs_HbO_zero_error = std_to_analyse(:,1);
    std_wavs_HbR_zero_error = std_to_analyse(:,2);
    std_wavs_CCO_zero_error = std_to_analyse(:,3);

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'SNR 5 wavs concs2.png', 'png')
shadedErrorBar(xvalues,change_wavs_CCO_zero_error,std_wavs_CCO_zero_error,'lineProps',{'g-o','markerfacecolor','g','markersize',10})
%shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
shadedErrorBar(xvalues,change_wavs_HbR_zero_error,std_wavs_HbR_zero_error,'lineProps',{'b-o','markerfacecolor','b','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbO_zero_error,std_wavs_HbO_zero_error,'lineProps',{'r-o','markerfacecolor','r','markersize',10})
grid on
ax = gca;
ax.FontSize = 40;
%xlim([0 200]);
ylim([0 1000]);
legend('oxCCO','HbR','HbO','Location','Best');
xlabel('Bandwidth');
ylabel({'Average percentage error in', 'chromophore concentration values'});
%title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
saveas(gcf,'SNR 5 wavs concs2.png')