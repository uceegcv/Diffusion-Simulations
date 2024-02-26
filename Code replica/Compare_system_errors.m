broad_errors = load('all_BB_errors.mat');

AV3_error = load('5wav_LED_error.mat');

laser_is_error = load('5wav_laser_Israel_error.mat');

percent_vector = [0.01; 5; 10; 15];

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing BB system errors HbO.png', 'png')
plot(percent_vector,broad_errors.ave_pert_error_BBSNR50(:,1),'rx-','MarkerSize',20,'LineWidth',1.8)
hold on
plot (percent_vector,broad_errors.ave_pert_error_BB2SNR50(:,1),'b-','MarkerSize',20,'LineWidth',1.8)
plot (percent_vector,broad_errors.ave_pert_error_BB3SNR50(:,1),'g-','MarkerSize',20,'LineWidth',1.8)
plot (percent_vector,broad_errors.ave_pert_error_BB4SNR50(:,1),'m-','MarkerSize',20,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
xlim([-2 17]);
legend('BB1','BB2','BB3','BB4','Location','Best');
xlabel('Percentage change in both extinction coefficient and DPF values');
ylabel('Average percentage error in chromophore concentration values');
title({'Comparing broadband system errors in HbO with varying error in DPF and E'});
saveas(gcf,'Comparing BB system errors HbO.png')
