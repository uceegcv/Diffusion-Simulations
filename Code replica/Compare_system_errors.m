broad_errors = load('all_BB_errors.mat');

AV3_error = load('5wav_LED_error.mat');

laser_is_error = load('5wav_laser_Israel_error.mat');

percent_vector = [0.01; 5; 10; 15];

% wavelengths_vector = [6 7 8 9 10 11 12 13 14 15 16 18 19 20 22 24 27 30 35 40 48 60 80 120 240];

xvalues = 1:10:201;

% filePattern = fullfile('change_wavelength_errors','*mat');

% theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/change_wavelengths_errors/*.mat');
%theFiles = dir('\\ifs.eng.cam.ac.uk\users\gcl33\Documents\GitHub\Diffusion-Simulations\Code replica\change_wavelengths_errors\*.mat');
% theFiles = dir('\\ifs.eng.cam.ac.uk\users\gcl33\Documents\GitHub\Diffusion-Simulations\Code replica\Bandwidth_files\*.mat');
theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Bandwidth find min noiseless original DPF/*.mat');

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
%     std_to_analyse = newdata.final_std_val_LEDSNR50;
    
    for d = 1:length(percent_vector)
        change_wavs_HbO_zero_error(k,d) = data_to_analyse(d,1);
        change_wavs_HbR_zero_error(k,d) = data_to_analyse(d,2);
        change_wavs_CCO_zero_error(k,d) = data_to_analyse(d,3);
%         std_wavs_HbO_zero_error(k,d) = std_to_analyse(d,1);
%         std_wavs_HbR_zero_error(k,d) = std_to_analyse(d,2);
%         std_wavs_CCO_zero_error(k,d) = std_to_analyse(d,3);
    end
end

%% Plot

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 1-200band systems zero parameter error shaded noiseless orig DPF.png', 'png')
plot(xvalues,change_wavs_HbO_zero_error(:,1),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
plot (xvalues,change_wavs_HbR_zero_error(:,1),'bx-','MarkerSize',10,'LineWidth',1.8)
plot (xvalues,change_wavs_CCO_zero_error(:,1),'gx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([20 50]);
legend('HbO','HbR','oxCCO','Location','Best');
xlabel('Bandwidth');
ylabel('Average percentage error in chromophore concentration values');
%title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
saveas(gcf,'Comparing 1-200band systems zero parameter error shaded noiseless orig DPF.png')

% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing BB systems zero parameter error.png', 'png')
% plot(wavelengths_vector,change_wavs_HbO_zero_error(:,1),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% plot (wavelengths_vector,change_wavs_HbR_zero_error(:,1),'bx-','MarkerSize',10,'LineWidth',1.8)
% plot (wavelengths_vector,change_wavs_CCO_zero_error(:,1),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([-2 17]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% title({'Comparing broadband system errors in chromophores with zero error in DPF and E'});
% saveas(gcf,'Comparing BB systems zero parameter error.png')
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
% imwrite(F.cdata, 'Comparing BB systems fifteen percent parameter error.png', 'png')
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

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 5wav 1-50BW systems zero parameter error no conv.png', 'png')
plot(xvalues,change_wavs_HbO_zero_error(:,1),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
plot (xvalues,change_wavs_HbR_zero_error(:,1),'bx-','MarkerSize',10,'LineWidth',1.8)
plot (xvalues,change_wavs_CCO_zero_error(:,1),'gx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([-2 17]);
legend('HbO','HbR','CCO','Location','Best');
xlabel('FWHM');
ylabel('Average percentage error in chromophore concentration values');
%title({'Comparing 5-wavelength systems with increasing bandwidth and zero error in DPF and E'});
saveas(gcf,'Comparing 5wav 1-50BW systems zero parameter error no conv.png')

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 5wav 1-50BW systems five percent parameter error no conv.png', 'png')
plot(xvalues,change_wavs_HbO_zero_error(:,2),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
plot (xvalues,change_wavs_HbR_zero_error(:,2),'bx-','MarkerSize',10,'LineWidth',1.8)
plot (xvalues,change_wavs_CCO_zero_error(:,2),'gx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([-2 17]);
legend('HbO','HbR','CCO','Location','Best');
xlabel('FWHM');
ylabel('Average percentage error in chromophore concentration values');
%title({'Comparing 5-wavelength systems with increasing bandwidth and five percent error in DPF and E'});
saveas(gcf,'Comparing 5wav 1-50BW systems five percent parameter error no conv.png')

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 5wav 1-50BW systems ten percent parameter error no conv.png', 'png')
plot(xvalues,change_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
plot (xvalues,change_wavs_HbR_zero_error(:,3),'bx-','MarkerSize',10,'LineWidth',1.8)
plot (xvalues,change_wavs_CCO_zero_error(:,3),'gx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([-2 17]);
legend('HbO','HbR','CCO','Location','Best');
xlabel('FWHM');
ylabel('Average percentage error in chromophore concentration values');
%title({'Comparing 5-wavelength systems with increasing bandwidth and ten percent error in DPF and E'});
saveas(gcf,'Comparing 5wav 1-50BW systems ten percent parameter error no conv.png')

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 5wav 1-50BW systems fifteen percent parameter error no conv.png', 'png')
plot(xvalues,change_wavs_HbO_zero_error(:,4),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
plot (xvalues,change_wavs_HbR_zero_error(:,4),'bx-','MarkerSize',10,'LineWidth',1.8)
plot (xvalues,change_wavs_CCO_zero_error(:,4),'gx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([-2 17]);
legend('HbO','HbR','CCO','Location','Best');
xlabel('FWHM');
ylabel('Average percentage error in chromophore concentration values');
%title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
saveas(gcf,'Comparing 5wav 1-50BW systems fifteen percent parameter error no conv.png')
