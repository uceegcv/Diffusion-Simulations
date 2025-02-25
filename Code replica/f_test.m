
f_test_data = xlsread('FandR2_240126.xlsx');

wavelength_column = f_test_data(:,1);

wavelength_vector = wavelength_column(9:end);

f_val_columns = f_test_data(:,13:16);

vals_to_plot = f_val_columns(9:end,:);


%%

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'f test four p vals.png', 'png')
plot(wavelength_vector(1:33),vals_to_plot(1:33,1),'x-','MarkerSize',10,'LineWidth',1.8)
hold on
plot(wavelength_vector(1:33),vals_to_plot(1:33,2),'x-','MarkerSize',10,'LineWidth',1.8)
plot(wavelength_vector(1:33),vals_to_plot(1:33,3),'x-','MarkerSize',10,'LineWidth',1.8)
plot(wavelength_vector(1:33),vals_to_plot(1:33,4),'x-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([5 10]);
legend('P = 0.10','P = 0.05','P = 0.025','P = 0.01','Location','Best');
xlabel('Number of wavelengths');
ylabel('Ratio of sum squared of residuals');
%title({'Comparing errors in chromophores with fifteen percent error in DPF and E'});
saveas(gcf,'f test four p vals.png')

%%
FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'f test four p vals 5-10 wavs.png', 'png')
plot(wavelength_vector(1:33),vals_to_plot(1:33,1),'x-','MarkerSize',20,'LineWidth',1.8)
hold on
plot(wavelength_vector(1:33),vals_to_plot(1:33,2),'x-','MarkerSize',20,'LineWidth',1.8)
plot(wavelength_vector(1:33),vals_to_plot(1:33,3),'x-','MarkerSize',20,'LineWidth',1.8)
plot(wavelength_vector(1:33),vals_to_plot(1:33,4),'x-','MarkerSize',20,'LineWidth',1.8)
grid on
ax = gca;
xlim([5 10]);
ax.FontSize = 20;
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
legend('P = 0.10','P = 0.05','P = 0.025','P = 0.01','Location','Best');
xlabel('Number of wavelengths');
ylabel('Ratio of sum squared of residuals');
%title({'Comparing errors in chromophores with fifteen percent error in DPF and E'});
saveas(gcf,'f test four p vals 5-10 wavs.png')


