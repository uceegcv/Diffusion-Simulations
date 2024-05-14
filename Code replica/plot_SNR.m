%% Setup initial values

% theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Cond_errors_3to50_best/*.mat');
% % theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Errors 1-200band with conv correct std/*.mat');
%  theFiles = dir('\\ifs.eng.cam.ac.uk\users\gcl33\Documents\GitHub\Diffusion-Simulations\Code replica\Errors 1-200band no conv concs2\*.mat');
% 
% percent_vector = [0.01; 5; 10; 15];
% percent_vector = 10;

xvalues = [20;30;40;50;60;100];

%% Sort the files in order of ascending number of wavelengths

% theFiles = theFiles(~[theFiles.isdir]);
% 
% for f = 1:length(theFiles)
%     filestosort{1,f} = theFiles(f).name;
% end
% 
% [sortedFiles,idx] = natsortfiles(filestosort);
% theFiles = theFiles(idx);
% 
% for k = 1 : length(theFiles)
%     baseFileName = theFiles(k).name;
    newdata = load('ave_error_5_3to200wavelengths_SNR50_origDPF_perfectE_SNRvals_concs2 .mat');
    data_to_analyse = newdata.ave_error_SNR;
    std_to_analyse = newdata.stand_dev_error_SNR;
%     data_to_analyse = newdata.ave_perfectE_error;
%     std_to_analyse = newdata.std_perfectE;
    %conc_to_analyse = newdata.concs2_saved;
    
%     for d = 1:length(percent_vector)
        change_wavs_HbO_zero_error = data_to_analyse(:,1);
        change_wavs_HbR_zero_error = data_to_analyse(:,2);
        change_wavs_CCO_zero_error = data_to_analyse(:,3);
        std_wavs_HbO_zero_error = std_to_analyse(:,1);
        std_wavs_HbR_zero_error = std_to_analyse(:,2);
        std_wavs_CCO_zero_error = std_to_analyse(:,3);
%     end

% for w = 1:length(conc_to_analyse)
%     HbO_conc(w) = conc_to_analyse{w}(1);
%     HbR_conc(w) = conc_to_analyse{w}(2);
%     CCO_conc(w) = conc_to_analyse{w}(3);
%     HbO_conc_change(w) = HbO_conc(w) - 56;
%     HbR_conc_change(w) = HbR_conc(w) - 24;
%     CCO_conc_change(w) = CCO_conc(w) - 4.9;
% end


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
%legend('oxCCO','HbR','HbO','Location','Best');
xlabel('Bandwidth');
ylabel({'Average percentage error in', 'chromophore concentration values'});
%title({'Comparing 5-wavelength systems with increasing bandwidth and fifteen percent error in DPF and E'});
saveas(gcf,'SNR 5 wavs concs2.png')