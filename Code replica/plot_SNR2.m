%% Setup initial values

% theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Cond_errors_3to50_best/*.mat');
 theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Errors SNR concs1/*.mat');
%  theFiles = dir('\\ifs.eng.cam.ac.uk\users\gcl33\Documents\GitHub\Diffusion-Simulations\Code replica\Errors SNR concs1\*.mat');

%percent_vector = [0.01; 5; 10; 15];
% percent_vector = 10;

xvalues = 3:200;
SNR_vals = [20;30;40;50;60;100];

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
    data_to_analyse = newdata.ave_error_SNR;
    std_to_analyse = newdata.stand_dev_error_SNR;
%     data_to_analyse = newdata.ave_perfectE_error;
%     std_to_analyse = newdata.std_perfectE;
    %conc_to_analyse = newdata.concs2_saved;
    
    for d = 1:length(SNR_vals)
        change_wavs_HbO_zero_error(k,d) = data_to_analyse(d,1);
        change_wavs_HbR_zero_error(k,d) = data_to_analyse(d,2);
        change_wavs_CCO_zero_error(k,d) = data_to_analyse(d,3);
        std_wavs_HbO_zero_error(k,d) = std_to_analyse(d,1);
        std_wavs_HbR_zero_error(k,d) = std_to_analyse(d,2);
        std_wavs_CCO_zero_error(k,d) = std_to_analyse(d,3);
    end
end


FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 3-200 wavs different SNRs concs1 CCO.png', 'png')
shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,1),std_wavs_CCO_zero_error(:,1),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
%shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,2),std_wavs_CCO_zero_error(:,2),'lineProps',{'g-o','markerfacecolor','g','markersize',10})
shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,3),std_wavs_CCO_zero_error(:,3),'lineProps',{'b-o','markerfacecolor','b','markersize',10})
shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,4),std_wavs_CCO_zero_error(:,4),'lineProps',{'m-o','markerfacecolor','m','markersize',10})
shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,5),std_wavs_CCO_zero_error(:,5),'lineProps',{'y-o','markerfacecolor',[0.8500 0.3250 0.0980],'markersize',10})
shadedErrorBar(xvalues,change_wavs_CCO_zero_error(:,6),std_wavs_CCO_zero_error(:,6),'lineProps',{'c-o','markerfacecolor','c','markersize',10})

%shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
grid on
ax = gca;
ax.FontSize = 38;
%xlim([0 200]);
ylim([0 7000]);
%legend('SNR 20','SNR 30','SNR 40','SNR 50','SNR 60','SNR 100','Location','Best');
xlabel('Number of Wavelengths');
ylabel({'Average percentage error in', 'chromophore concentration values'});
title({'oxCCO'});
saveas(gcf,'Comparing 3-200 wavs different SNRs concs1 CCO.png')


FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 3-200 wavs different SNRs concs1 HbO.png', 'png')
shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,1),std_wavs_HbO_zero_error(:,1),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
%shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,2),std_wavs_HbO_zero_error(:,2),'lineProps',{'g-o','markerfacecolor','g','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'lineProps',{'b-o','markerfacecolor','b','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,4),std_wavs_HbO_zero_error(:,4),'lineProps',{'m-o','markerfacecolor','m','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,5),std_wavs_HbO_zero_error(:,5),'lineProps',{'y-o','markerfacecolor',[0.8500 0.3250 0.0980],'markersize',10})
shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,6),std_wavs_HbO_zero_error(:,6),'lineProps',{'c-o','markerfacecolor','c','markersize',10})

%shadedErrorBar(xvalues,change_wavs_HbO_zero_error(:,3),std_wavs_HbO_zero_error(:,3),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
grid on
ax = gca;
ax.FontSize = 38;
%xlim([0 200]);
ylim([0 5000]);
legend('SNR 20','SNR 30','SNR 40','SNR 50','SNR 60','SNR 100','Location','Best');
xlabel('Number of Wavelengths');
ylabel({'Average percentage error in', 'chromophore concentration values'});
title({'HbO'});
saveas(gcf,'Comparing 3-200 wavs different SNRs concs1 HbO.png')

FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparing 3-200 wavs different SNRs concs1 HbR.png', 'png')
shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,1),std_wavs_HbR_zero_error(:,1),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
%shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,3),std_wavs_HbR_zero_error(:,3),'rx-','MarkerSize',10,'LineWidth',1.8)
hold on
shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,2),std_wavs_HbR_zero_error(:,2),'lineProps',{'g-o','markerfacecolor','g','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,3),std_wavs_HbR_zero_error(:,3),'lineProps',{'b-o','markerfacecolor','b','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,4),std_wavs_HbR_zero_error(:,4),'lineProps',{'m-o','markerfacecolor','m','markersize',10})
shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,5),std_wavs_HbR_zero_error(:,5),'lineProps',{'y-o','markerfacecolor',[0.8500 0.3250 0.0980],'markersize',10})
shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,6),std_wavs_HbR_zero_error(:,6),'lineProps',{'c-o','markerfacecolor','c','markersize',10})

%shadedErrorBar(xvalues,change_wavs_HbR_zero_error(:,3),std_wavs_HbR_zero_error(:,3),'lineProps',{'r-o','markerfacecolor','r','markersize',10})
grid on
ax = gca;
ax.FontSize = 38;
%xlim([0 200]);
ylim([0 3000]);
%legend('SNR 20','SNR 30','SNR 40','SNR 50','SNR 60','SNR 100','Location','Best');
xlabel('Number of Wavelengths');
ylabel({'Average percentage error in', 'chromophore concentration values'});
title({'HbR'});
saveas(gcf,'Comparing 3-200 wavs different SNRs concs1 HbR.png')
