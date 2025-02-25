%% Setup initial values

% theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Cond_errors_3to50_best/*.mat');
 theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Errors SNR concs1/*.mat');
% theFiles = dir('\\ifs.eng.cam.ac.uk\users\gcl33\Documents\GitHub\Diffusion-Simulations\Code replica\Errors 3-200wav correct std\*.mat');

percent_vector = [0.01; 5; 10; 15];
% percent_vector = 10;

xvalues = 3:200;

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
%     data_to_analyse = newdata.ave_perfectE_error;
%     std_to_analyse = newdata.std_perfectE;
    %conc_to_analyse = newdata.concs2_saved;
    
    % For when there are four values for different % uncertainties
    for d = 1:length(percent_vector)
        change_wavs_HbO_zero_error(k,d) = data_to_analyse(d,1);
        change_wavs_HbR_zero_error(k,d) = data_to_analyse(d,2);
        change_wavs_CCO_zero_error(k,d) = data_to_analyse(d,3);
        std_wavs_HbO_zero_error(k,d) = std_to_analyse(d,1);
        std_wavs_HbR_zero_error(k,d) = std_to_analyse(d,2);
        std_wavs_CCO_zero_error(k,d) = std_to_analyse(d,3);
    end

%         % For perfect E dataset
%         change_wavs_HbO_zero_error(k,1) = data_to_analyse(1,1);
%         change_wavs_HbR_zero_error(k,1) = data_to_analyse(1,2);
%         change_wavs_CCO_zero_error(k,1) = data_to_analyse(1,3);
%         std_wavs_HbO_zero_error(k,1) = std_to_analyse(1,1);
%         std_wavs_HbR_zero_error(k,1) = std_to_analyse(1,2);
%         std_wavs_CCO_zero_error(k,1) = std_to_analyse(1,3);
end

% for w = 1:length(conc_to_analyse)
%     HbO_conc(w) = conc_to_analyse{w}(1);
%     HbR_conc(w) = conc_to_analyse{w}(2);
%     CCO_conc(w) = conc_to_analyse{w}(3);
%     HbO_conc_change(w) = HbO_conc(w) - 56;
%     HbR_conc_change(w) = HbR_conc(w) - 24;
%     CCO_conc_change(w) = CCO_conc(w) - 4.9;
% end

% % For regular wav graphs
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing best combos 3-50wav systems ten percent parameter error.png', 'png')
% plot(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% errorbar(xvalues,change_wavs_HbR_zero_error(:,3),std_wavs_HbR_zero_error(:,3),std_wavs_HbR_zero_error(:,3),'bx-','MarkerSize',10,'LineWidth',1.8)
% errorbar(xvalues,change_wavs_CCO_zero_error(:,3),std_wavs_CCO_zero_error(:,3),std_wavs_CCO_zero_error(:,3),'gx-','MarkerSize',10,'LineWidth',1.8)
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([20 50]);
% legend('HbO','HbR','CCO','Location','Best');
% xlabel('Number of wavelengths');
% ylabel('Average percentage error in chromophore concentration values');
% %title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing best combos 3-50wav systems ten percent parameter error.png')

% % For perfect E
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Comparing 3-200 wavs different Es paper.png', 'png')
% shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,1),std_wavs_CCO_zero_error(:,1),'lineProps',{'g-o','markerfacecolor','g','markersize',10})
% %shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
% hold on
% shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,1),std_wavs_HbR_zero_error(:,1),'lineProps',{'b-o','markerfacecolor','b','markersize',10})
% shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,1),std_wavs_HbO_zero_error(:,1),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
% grid on
% ax = gca;
% ax.FontSize = 20;
% %xlim([0 200]);
% %ylim([0 35]);
% legend('oxCCO','HbR','HbO','Location','Best');
% xlabel('Number of Wavelengths');
% % ylabel({'Average percentage error in', 'chromophore concentration values'});
% ylabel({'Average percentage error in chromophore concentration values'});
% 
% %title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
% saveas(gcf,'Comparing 3-200 wavs different Es paper.png')


% For bandwidth graphs
FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 3-200wav errors concs1 paper.png', 'png')
shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,3),std_wavs_CCO_zero_error(:,3),'lineProps',{'g-o','markerfacecolor','g','markersize',10})
%shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,3),std_wavs_HbR_zero_error(:,3),'lineProps',{'b-o','markerfacecolor','b','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
grid on
ax = gca;
ax.FontSize = 20;
%xlim([0 70]);
ylim([0 600]);
%ylim([0 50]);
legend('oxCCO','HbR','HbO','Location','Best');
xlabel('Number of Wavelengths');
% ylabel({'Average percentage error in', 'chromophore concentration values'});
ylabel('Average percentage error in chromophore concentration values');

%title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
saveas(gcf,'Comparing 3-200wav errors concs1 paper.png')

%% Correlation plot

for p = 1:length(change_wavs_HbO_zero_error(:,1))
    correlation(p) = (change_wavs_HbO_zero_error(p,1) + change_wavs_CCO_zero_error(p,1))/2;
end

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 3-200 wav systems zero parameter error zeroed correlation.png', 'png')
plot(xvalues, correlation,'rx-','MarkerSize',10,'LineWidth',1.8)
grid on
ax = gca;
ax.FontSize = 20;
%xlim([0 50]);
%ylim([0 10]);
%legend('oxCCO','HbR','HbO','Location','Best');
xlabel('Number of wavelengths');
ylabel('Correlation');
%title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
saveas(gcf,'Comparing 3-200 wav systems zero parameter error zeroed correlation.png')



