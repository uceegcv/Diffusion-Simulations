broad_errors = load('all_BB_errors.mat');

% AV3_error = load('5wav_LED_error.mat');
AV3_error = load('AV3_error_wavelengths_std.mat');

laser_is_error = load('5wav_laser_Israel_error.mat');

percent_vector = [0.01; 5; 10; 15];

% wavelengths_vector = [5 6 7 8 9 10 11 12 13 14 15 16 18 19 20 22 24 27 30 35 40 48 60 80 120 240];
wavelengths_vector = 5:200;

filePattern = fullfile('change_wavelength_errors','*mat');

% theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Ave_std_files/*.mat');
%theFiles = dir('\\ifs.eng.cam.ac.uk\users\gcl33\Documents\GitHub\Diffusion-Simulations\Code replica\change_wavelengths_errors\*.mat');
theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Full wavelength range/*.mat');

%% Sort the files in order of ascending number of wavelengths

theFiles = theFiles(~[theFiles.isdir]);

for f = 1:length(theFiles)
    filestosort{1,f} = theFiles(f).name;
end

[sortedFiles,idx] = natsortfiles(filestosort);
theFiles = theFiles(idx);

for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    newdata = load(baseFileName);
    data_to_analyse = newdata.ave_pert_error_LEDSNR50;
    std_to_analyse = newdata.final_std_val_LEDSNR50;
    
    for d = 1:length(percent_vector)
        change_wavs_HbO(k,d) = data_to_analyse(d,1);
        change_wavs_HbR(k,d) = data_to_analyse(d,2);
        change_wavs_CCO(k,d) = data_to_analyse(d,3);
        std_wavs_HbO(k,d) = std_to_analyse(d,1);
        std_wavs_HbR(k,d) = std_to_analyse(d,2);
        std_wavs_CCO(k,d) = std_to_analyse(d,3);
    end
end

% Add in the five-wavelength error from the AV3 file
% for r = 1:4
%     five_wav_HbO(r) = AV3_error.ave_pert_error_LEDSNR50(r,1);
%     five_wav_HbR(r) = AV3_error.ave_pert_error_LEDSNR50(r,2);
%     five_wav_CCO(r) = AV3_error.ave_pert_error_LEDSNR50(r,3);
%     five_wav_HbO_std(r) = AV3_error.final_std_val_LEDSNR50(r,1);
%     five_wav_HbR_std(r) = AV3_error.final_std_val_LEDSNR50(r,2);
%     five_wav_CCO_std(r) = AV3_error.final_std_val_LEDSNR50(r,3);
% end
%     
% change_wavs_HbO_zero_error = [five_wav_HbO;change_wavs_HbO];
% change_wavs_HbR_zero_error = [five_wav_HbR;change_wavs_HbR];
% change_wavs_CCO_zero_error = [five_wav_CCO;change_wavs_CCO];
% change_wavs_HbO_zero_std = [five_wav_HbO_std;std_wavs_HbO];
% change_wavs_HbR_zero_std = [five_wav_HbR_std;std_wavs_HbR];
% change_wavs_CCO_zero_std = [five_wav_CCO_std;std_wavs_CCO];

%% Plot
% 
FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 5-200 wavs zero parameter error.png', 'png')
errorbar(wavelengths_vector,change_wavs_HbO(:,1),std_wavs_HbO(:,1),std_wavs_HbO(:,1),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
errorbar(wavelengths_vector,change_wavs_HbR(:,1),std_wavs_HbR(:,1),std_wavs_HbR(:,1),'bx-','MarkerSize',10,'LineWidth',1.8)
errorbar(wavelengths_vector,change_wavs_CCO(:,1),std_wavs_CCO(:,1),std_wavs_CCO(:,1),'gx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([-2 17]);
legend('HbO','HbR','CCO','Location','Best');
xlabel('Number of wavelengths');
ylabel('Average percentage error in chromophore concentration values');
%title({'Comparing broadband system errors in chromophores with zero error in DPF and E'});
saveas(gcf,'Comparing 5-200 wavs zero parameter error.png')

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 5-200 wavs ten percent parameter error.png', 'png')
errorbar(wavelengths_vector,change_wavs_HbO(:,3),std_wavs_HbO(:,3),std_wavs_HbO(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
errorbar(wavelengths_vector,change_wavs_HbR(:,3),std_wavs_HbR(:,3),std_wavs_HbR(:,3),'bx-','MarkerSize',10,'LineWidth',1.8)
errorbar(wavelengths_vector,change_wavs_CCO(:,3),std_wavs_CCO(:,3),std_wavs_CCO(:,3),'gx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([-2 17]);
legend('HbO','HbR','CCO','Location','Best');
xlabel('Number of wavelengths');
ylabel('Average percentage error in chromophore concentration values');
%title({'Comparing broadband system errors in chromophores with zero error in DPF and E'});
saveas(gcf,'Comparing 5-200 wavs ten percent parameter error.png')
% 
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing BB systems five percent parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,2),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,2),'bx-','MarkerSize',10,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,2),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([-2 17]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% title({'Comparing broadband system errors in chromophores with five percent error in DPF and E'});
% saveas(gcf,'Comparing BB systems five percent parameter error.png')
% 
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing BB systems ten percent parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,3),'bx-','MarkerSize',10,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,3),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([-2 17]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% title({'Comparing broadband system errors in chromophores with ten percent error in DPF and E'});
% saveas(gcf,'Comparing BB systems ten percent parameter error.png')
% 
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing BB systems fiteen percent parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,4),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,4),'bx-','MarkerSize',10,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,4),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([-2 17]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% title({'Comparing broadband system errors in chromophores with fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing BB systems fifteen percent parameter error.png')
% 
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing wavs 5-10 fiteen percent parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,4),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,4),'bx-','MarkerSize',10,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,4),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% xlim([5 10]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% %title({'Comparing errors in chromophores with fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing wavs 5-10 fiteen percent parameter error.png')
% %%
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing wavs 5-8 fiteen percent parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,4),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,4),'bx-','MarkerSize',10,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,4),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% xlim([5 8]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% %title({'Comparing errors in chromophores with fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing wavs 5-8 fiteen percent parameter error.png')
% 
% %%
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing wavs 5-10 fifteen percent parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,4),'rx-','MarkerSize',20,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,4),'bx-','MarkerSize',20,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,4),'gx-','MarkerSize',20,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% xlim([5 10]);
% curtick = get(gca, 'xTick');
% xticks(unique(round(curtick)));
% %ax.xTicks = unique(round(ax.xTicks));
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% %title({'Comparing errors in chromophores with fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing wavs 5-10 fifteen percent parameter error.png')
% 
% %%
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing wavs 5-242 fifteen percent parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,4),'rx-','MarkerSize',20,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,4),'bx-','MarkerSize',20,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,4),'gx-','MarkerSize',20,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([5 10]);
% curtick = get(gca, 'xTick');
% xticks(unique(round(curtick)));
% %ax.xTicks = unique(round(ax.xTicks));
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% %title({'Comparing errors in chromophores with fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing wavs 5-242 fifteen percent parameter error.png')

% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing all wavs fiteen percent parameter error.png', 'png')
% errorbar(wavelengths_vector,change_wavs_HbO_zero_error(:,4),change_wavs_HbO_zero_std(:,4),change_wavs_HbO_zero_std(:,4),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% errorbar(wavelengths_vector,change_wavs_HbR_zero_error(:,4),change_wavs_HbR_zero_std(:,4),change_wavs_HbR_zero_std(:,4),'bx-','MarkerSize',10,'LineWidth',1.8)
% errorbar(wavelengths_vector,change_wavs_CCO_zero_error(:,4),change_wavs_CCO_zero_std(:,4),change_wavs_CCO_zero_std(:,4),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([5 10]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% %title({'Comparing errors in chromophores with fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing all wavs fiteen percent parameter error.png')

