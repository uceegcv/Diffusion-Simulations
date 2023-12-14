%% Diffuse reflection model for a slab geometry
% % September 2023
% clear all;
% clc;
% close all;
% addpath(genpath('N:\Matlab code\homer2_v2_8_11022018')) 

% Enter simulation data
SD = 30; % Source-detector separation on same surface of slab (mm)
slab = 50; % Slab thickness (mm)

%% Code for Broadband System 1
% x contains all broadband wavelengths
x = 740:921;

% Contains broadband wavelengths of interest - the one to edit
%wl = 680:921; %nm (cannot be a larger range than x)
wl = 740:921; % de Roever

% Load chromophore data
% Generate vectors containing extinction coefficient and mua data
% Using Journal of Biomedical Optics 15 5, 056002 September/October 2010
% GetExtinctions(lambda)
% Returns the specific absorption coefficients for [HbO Hb H2O lipid aa3]
% for the specified wavelengths. Note that the specific
% absorption coefficient (defined base e) is equal to the 
% specific extinction coefficient (defined base 10) times 2.303.
E = GetExtinctions(x);
E_new = GetExtinctions(wl);

%Yields specific absorption coefficients in units of cm-1 per Molar or 
%absorption coefficient in cm-1 for water
%Convert to mm-1 per uM.
%1 cm-1 per Molar is 1e-1 mm-1 per Molar, which is 1e-7 mm-1 per microMolar
E(:,1) = E(:,1)/1e7; %E_HbO
E(:,2) = E(:,2)/1e7; %E_HbR
E(:,3) = E(:,3)/10; %mua_H2O
E(:,4) = E(:,4)/10; %mua_lipid
E(:,5) = E(:,5)/1e7; %E_aa3

E_new(:,1) = E_new(:,1)/1e7; %E_HbO
E_new(:,2) = E_new(:,2)/1e7; %E_HbR
E_new(:,3) = E_new(:,3)/10; %mua_H2O
E_new(:,4) = E_new(:,4)/10; %mua_lipid
E_new(:,5) = E_new(:,5)/1e7; %E_aa3

%Replace lipid values from those on OMLC database (and convert to mm-1)
Lvalue = Lipid_mua(x)/1000;
E(:,4) = Lvalue(:);
Lvalue2 = Lipid_mua(wl)/1000;
E_new(:,4) = Lvalue2(:);

E_rewrite = [E(:,1) E(:,2) E(:,5)]; % ADD
E = E_rewrite;

E_new_rewrite = [E_new(:,1) E_new(:,2) E_new(:,5)]; % ADD
E_new = E_new_rewrite;

%% Code for Broadband System

W(1) = 0.8;    %Water fraction
L1(1) = 0.116;  %Lipid fraction
B(1) = 0.012;  %Background non-wavelength-dependent absorption coefficient (mm-1)
C_HbO(1) = 56; %Concentration of HbO (uM)
C_HbR(1) = 24; %Concentration of HbR (uM)
C_aa3(1) = 4.9; %Concentration of aa3 (uM)
% mua1 = (E(:,1)*C_HbO(1)) + (E(:,2)*C_HbR(1)) + (E(:,3)*W(1)) + (E(:,4)*L1(1)) + (E(:,5)*C_aa3(1));
mua1 = (E(:,1)*C_HbO(1)) + (E(:,2)*C_HbR(1)) + (E(:,3)*C_aa3(1));

% Put concs in vector
% concs1 = [C_HbO(1) C_HbR(1) W(1) L1(1) C_aa3(1)];
concs1 = [C_HbO(1) C_HbR(1) C_aa3(1)];

I_original = exp(-mua1);

%% Change concentrations
    
W(2) = 0.8;    %Water fraction
L1(2) = 0.116;  %Lipid fraction
B(2) = 0.012;  %Background non-wavelength-dependent absorption coefficient (mm-1)
C_HbO(2) = C_HbO(1) + 1.5; %Concentration of HbO (uM)
C_HbR(2) = C_HbR(1) - 0.5; %Concentration of HbR (uM)
C_aa3(2) = C_aa3(1) + 0.4; %Concentration of aa3 (uM)
% mua2 = (E(:,1)*C_HbO(2)) + (E(:,2)*C_HbR(2)) + (E(:,3)*W(2)) + (E(:,4)*L1(2)) + (E(:,5)*C_aa3(2));
mua2 = (E(:,1)*C_HbO(2)) + (E(:,2)*C_HbR(2)) + (E(:,3)*C_aa3(2));

% concs2 = [C_HbO(2) C_HbR(2) W(2) L1(2) C_aa3(2)];
concs2 = [C_HbO(2) C_HbR(2) C_aa3(2)];

% Define actual (hard-coded) change in concs
ac_conc_change = concs2 - concs1;
ac_conc_change = ac_conc_change.';
    
%% Generate a typical brain spectrum

    % Typical chromophone concentration

%     figure;
%     plot(x,mua1,'ko-');
%     hold on
%     plot(x,mua2,'ko-');
%     xlim([min(x) max(x)])
%     xlabel('Wavelength (nm)');
%     ylabel('Absorption coefficient mm^{-1})');
%     title('Tissue absorption spectrum');
% 
%     % Define and plot transport scatter coefficient
    power = -1.2;
    musp = 1.0 * (x/800).^power;
%     figure
%     plot(x,musp,'ko-');
%     xlim([min(x) max(x)])
%     xlabel('Wavelength (nm)');
%     ylabel('Transport scatter coefficient (mm^{-1})');
%     title('Transport scatter coefficient');

    %% Slab model from Contini et al. Appl.Opt. 36, 4587 (1997). - 1
    
    mut = mua1 + musp';
    D = 1 ./ (3 * mut');
    mueff = sqrt(3 .* mua1' .* mut');
    z0 = 1 ./ musp;

    % Mismatched boundary parameters
    nrel = 1.4; % Refractive index of medium
    Apar = 504.332889-(2641.00214*nrel)+(5923.699064*nrel^2)-(7376.335814*nrel^3)+(5507.53041*nrel^4)-(2463.357945*nrel^5)+(610.956547*nrel^6)-(64.8047*nrel^7);
    ze = (2 * Apar) .* D;

    % Sum terms to generate reflectance: Equation (45)
    mlimit = 4;
    R1_BB = zeros([1 length(x)]); 
    %R(l) = num2cell(R);
    for m =(-mlimit):1:(mlimit)
        z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
        z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
        g = (SD * SD) + (z3m .* z3m);
        h = (SD * SD) + (z3m .* z3m);
        rmg = sqrt(mua1' .* g ./ D);
        rmh = sqrt(mua1' .* h ./ D);
        ermg = exp(-rmg);
        ermh = exp(-rmh);
        RT11 = z3m .* (g.^-1.5) .* (rmg + 1) .* ermg;
        RT21 = z4m .* (h.^-1.5) .* (rmh + 1) .* ermh;
        R1_BB = R1_BB + ((RT11 - RT21) ./ (-4 * pi));
%         R{l} = R;
    end
        %R{l} = num2cell(R);
    %R1(1,:) = R1;
    % Sum terms to generate DPF: Equation (47)
    mlimit = 4;
    L1_BB = zeros([1 length(x)]); 
    for m =(-mlimit):1:(mlimit)
        z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
        z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
        g = (SD * SD) + (z3m .* z3m);
        h = (SD * SD) + (z3m .* z3m);
        rmg = sqrt(mua1' .* g ./ D);
        rmh = sqrt(mua1' .* h ./ D);
        ermg = exp(-rmg);
        ermh = exp(-rmh);
        LT0 = (-8 * pi) .* D .* R1_BB;
        LT11 = z3m .* (g.^-0.5) .* ermg;
        LT21 = z4m .* (h.^-0.5) .* ermh;
        L1_BB = L1_BB + ((LT11 - LT21) ./ LT0); 
    end
    

%     %% Calculate and plot the DPF
    DPF1_BB = L1_BB ./ SD; % Wavelength dependent
%     figure
%     plot(x,DPF1_BB,'ko-');
%     xlim([min(x) max(x)])
%     xlabel('Wavelength (nm)');
%     ylabel('DPF');
%     title('DPF 1');
% 
%     %% Plot attenuation (-log diffuse reflectance)
    A1 = -log(R1_BB);
%     figure
%     subplot(1,2,1);
%     plot(x,R1_BB,'bo-');
%     xlim([min(x) max(x)]);
%     xlabel('Wavelength (nm)');
%     ylabel('Diffuse Reflectance');
%     title('Diffuse Reflectance 1');
%     subplot(1,2,2);
%     plot(x,A1,'ko-');
%     xlim([min(x) max(x)]);
%     xlabel('Wavelength (nm)');
%     ylabel('-log_{e}Reflectance');
%     title('Attenuation = -ln(Reflectance) 1');

    %% Slab model from Contini et al. Appl.Opt. 36, 4587 (1997). - 2
    
    mut = mua2 + musp';
    D = 1 ./ (3 * mut');
    mueff = sqrt(3 .* mua2' .* mut');
    z0 = 1 ./ musp;

    % Mismatched boundary parameters
    nrel = 1.4; % Refractive index of medium
    Apar = 504.332889-(2641.00214*nrel)+(5923.699064*nrel^2)-(7376.335814*nrel^3)+(5507.53041*nrel^4)-(2463.357945*nrel^5)+(610.956547*nrel^6)-(64.8047*nrel^7);
    ze = (2 * Apar) .* D;

    % Sum terms to generate reflectance: Equation (45)
    mlimit = 4;
    R2_BB = zeros([1 length(x)]); 
    %R(l) = num2cell(R);
    for m =(-mlimit):1:(mlimit)
        z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
        z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
        g = (SD * SD) + (z3m .* z3m);
        h = (SD * SD) + (z3m .* z3m);
        rmg = sqrt(mua2' .* g ./ D);
        rmh = sqrt(mua2' .* h ./ D);
        ermg = exp(-rmg);
        ermh = exp(-rmh);
        RT12 = z3m .* (g.^-1.5) .* (rmg + 1) .* ermg;
        RT22 = z4m .* (h.^-1.5) .* (rmh + 1) .* ermh;
        R2_BB = R2_BB + ((RT12 - RT22) ./ (-4 * pi));
%         R{l} = R;
    end
        %R{l} = num2cell(R);
    % Sum terms to generate DPF: Equation (47)
    mlimit = 4;
    L2_BB = zeros([1 length(x)]); 
    for m =(-mlimit):1:(mlimit)
        z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
        z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
        g = (SD * SD) + (z3m .* z3m);
        h = (SD * SD) + (z3m .* z3m);
        rmg = sqrt(mua2' .* g ./ D);
        rmh = sqrt(mua2' .* h ./ D);
        ermg = exp(-rmg);
        ermh = exp(-rmh);
        LT0 = (-8 * pi) .* D .* R2_BB;
        LT12 = z3m .* (g.^-0.5) .* ermg;
        LT22 = z4m .* (h.^-0.5) .* ermh;
        L2_BB = L2_BB + ((LT12 - LT22) ./ LT0); 
    end
    

%     %% Calculate and plot the DPF
    DPF2_BB = L2_BB ./ SD; % Wavelength dependent
%     figure
%     plot(x,DPF2_BB,'ko-');
%     xlim([min(x) max(x)])
%     xlabel('Wavelength (nm)');
%     ylabel('DPF');
%     title('DPF 2');
% 
%     %% Plot attenuation (-log diffuse reflectance)
    A2 = -log(R2_BB);
%     figure
%     subplot(1,2,1);
%     plot(x,R2_BB,'bo-');
%     xlim([min(x) max(x)]);
%     xlabel('Wavelength (nm)');
%     ylabel('Diffuse Reflectance');
%     title('Diffuse Reflectance 2');
%     subplot(1,2,2);
%     plot(x,A2,'ko-');
%     xlim([min(x) max(x)]);
%     xlabel('Wavelength (nm)');
%     ylabel('-log_{e}Reflectance');
%     title('Attenuation = -ln(Reflectance) 2');

%% Task 2

change_I_BB_diff = log(R2_BB) - log(R1_BB);

ave_DPF_BB_diff = (DPF1_BB + DPF2_BB)/2;

for t = 1:length(x)
    calc_I_BB(t) = ave_DPF_BB_diff(t)*SD*(C_HbO(2)-C_HbO(1))*E(t,1);
end

for t = 1:length(x)
    perc_diff_I(t) = 100 - ((calc_I_BB(t)/change_I_BB_diff(t))*100);
end

E_prime_invert_BB = pinv(E);

% This is the value to pay attention to!! From diffusion
conc_BB_diff =  E_prime_invert_BB * (-change_I_BB_diff ./ (SD*ave_DPF_BB_diff))'; % From diffusion

%%
for t = 1:length(ac_conc_change)
    perc_diff_BB_diff(t) = ((conc_BB_diff(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
end

SNR_vector = [20; 30; 40; 50; 60];

for l = 1:20 % Number of iterations to average for error bars

    R1_BB_SNR20 = awgn(R1_BB,20,'measured');
    R2_BB_SNR20 = awgn(R2_BB,20,'measured');
    change_I_BB_SNR20 = log(R2_BB_SNR20) - log(R1_BB_SNR20);
    % change_I_BB_SNR20 = log(R2_BB_SNR20) - log(R1_BB);

    R1_BB_SNR30 = awgn(R1_BB,30,'measured');
    R2_BB_SNR30 = awgn(R2_BB,30,'measured');
    change_I_BB_SNR30 = log(R2_BB_SNR30) - log(R1_BB_SNR30);
    % change_I_BB_SNR30 = log(R2_BB_SNR30) - log(R1_BB);

    R1_BB_SNR40 = awgn(R1_BB,40,'measured');
    R2_BB_SNR40 = awgn(R2_BB,40,'measured');
    change_I_BB_SNR40 = log(R2_BB_SNR40) - log(R1_BB_SNR40);
    % change_I_BB_SNR40 = log(R2_BB_SNR40) - log(R1_BB);

    R1_BB_SNR50 = awgn(R1_BB,50,'measured');
    R2_BB_SNR50 = awgn(R2_BB,50,'measured');
    change_I_BB_SNR50 = log(R2_BB_SNR50) - log(R1_BB_SNR50);
    % change_I_BB_SNR50 = log(R2_BB_SNR50) - log(R1_BB);

    R1_BB_SNR60 = awgn(R1_BB,60,'measured');
    R2_BB_SNR60 = awgn(R2_BB,60,'measured');
    change_I_BB_SNR60 = log(R2_BB_SNR60) - log(R1_BB_SNR60);
    % change_I_BB_SNR60 = log(R2_BB_SNR60) - log(R1_BB);

%     conc_BB_SNR20 =  E_prime_invert_BB * (-change_I_BB_SNR20 ./ (SD*ave_DPF_BB_diff))';
%     conc_BB_SNR30 =  E_prime_invert_BB * (-change_I_BB_SNR30 ./ (SD*ave_DPF_BB_diff))';
%     conc_BB_SNR40 =  E_prime_invert_BB * (-change_I_BB_SNR40 ./ (SD*ave_DPF_BB_diff))';
%     conc_BB_SNR50 =  E_prime_invert_BB * (-change_I_BB_SNR50 ./ (SD*ave_DPF_BB_diff))';
%     conc_BB_SNR60 =  E_prime_invert_BB * (-change_I_BB_SNR60 ./ (SD*ave_DPF_BB_diff))';
% 
% 
%     for t = 1:length(ac_conc_change)
%         perc_diff_BB(t) = ((conc_BB_diff(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_diff_BB_SNR20(t) = ((conc_BB_SNR20(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_diff_BB_SNR30(t) = ((conc_BB_SNR30(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_diff_BB_SNR40(t) = ((conc_BB_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_diff_BB_SNR50(t) = ((conc_BB_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_diff_BB_SNR60(t) = ((conc_BB_SNR60(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%     end
% 
%     perc_vector_SNR = [20;30;40;50;60;100];
%     perc_diff_SNR_vector = [perc_diff_BB_SNR20;perc_diff_BB_SNR30;perc_diff_BB_SNR40;perc_diff_BB_SNR50;perc_diff_BB_SNR60;perc_diff_BB];
%     perc_diff_SNR_vector = abs(perc_diff_SNR_vector);

%     perc_diff_SNR_vector_saved{:,l} = perc_diff_SNR_vector(:,:); % Save each iteration to a cell

    % FigH = figure('Position', get(0, 'Screensize'));
    % F = getframe(FigH);
    % imwrite(F.cdata, 'Comparison of percentage errors for increasing SNR of broadband system 1.png', 'png')
    % plot(perc_vector_SNR,perc_diff_SNR_vector(:,1),'r*-','MarkerSize',20)
    % hold on
    % plot(perc_vector_SNR,perc_diff_SNR_vector(:,2),'bs-','MarkerSize',10)
    % plot(perc_vector_SNR,perc_diff_SNR_vector(:,5),'go-','MarkerSize',10)
    % ax = gca;
    % ax.FontSize = 20;
    % legend('HbO','HbR','oxCCO','Location','Best');
    % xlabel('Absolute SNR value');
    % ylabel('Percentage error in chromophore concentration values');
    % title({'Comparison of percentage errors for increasing SNR', 'of broadband system 1 (680-921nm)'});
    % saveas(gcf,'Comparison of percentage errors for increasing SNR of broadband system 1.png')

    %% Introduce different DPF datasets
    
    % CYRIL dataset
    DPF_dep_cyril = stretch_DPF_Lambda_Dependency_680to915; 

    for r = 1:length(wl)
        DPF_index_BB(r) = find(DPF_dep_cyril(:,1) == wl(r).');
        DPF_vals_BB(r) = DPF_dep_cyril(DPF_index_BB(r),2);
    end

%     set_DPF = 4.99; % Infant
    set_DPF_cyril = 6; % Adult

    DPF_assumed_cyril = DPF_vals_BB * set_DPF_cyril;

    conc_cyril_DPF_BB = E_prime_invert_BB * (-change_I_BB_diff ./ (SD*DPF_assumed_cyril))';
    
    for t = 1:length(ac_conc_change)
        perc_diff_BB_cyrilDPF(t) = ((conc_cyril_DPF_BB(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
    end
    
    % Toolboxes with Scholkmann DPF vals
    
    alpha = 223.3;
    beta = 0.05624;
    gamma = 0.8493;
    sigma = -5.723e-7;
    epsilon = 0.001245;
    zeta = -0.9025;
    age1 = 20;
    age2 = 40;

    for i = 1:length(wl)
        DPF_scholk_age1(i) = alpha + (beta*(age1^gamma)) + (sigma*(wl(i)^3)) + (epsilon*(wl(i)^2)) + (zeta*wl(i));
    end
    
    for i = 1:length(wl)
        DPF_scholk_age2(i) = alpha + (beta*(age2^gamma)) + (sigma*(wl(i)^3)) + (epsilon*(wl(i)^2)) + (zeta*wl(i));
    end
    
    conc_scholk1_DPF_BB = E_prime_invert_BB * (-change_I_BB_diff ./ (SD*DPF_scholk_age1))';
    conc_scholk2_DPF_BB = E_prime_invert_BB * (-change_I_BB_diff ./ (SD*DPF_scholk_age2))';

    for t = 1:length(ac_conc_change)
        perc_diff_BB_scholk1(t) = ((conc_scholk1_DPF_BB(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
    end
    
    for t = 1:length(ac_conc_change)
        perc_diff_BB_scholk2(t) = ((conc_scholk2_DPF_BB(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
    end
   
    % DPF from homer2 - user input!!
    
    for i = 1:length(wl)
        DPF_homer(i) = 6;
    end
    
    conc_homer_DPF_BB = E_prime_invert_BB * (-change_I_BB_diff ./ (SD*DPF_homer))';
    
    for t = 1:length(ac_conc_change)
        perc_diff_BB_homer(t) = ((conc_homer_DPF_BB(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
    end
    
    FigH = figure('Position', get(0, 'Screensize'));
    F = getframe(FigH);
    imwrite(F.cdata, 'Comparison of DPF values used in different toolboxes.png', 'png')
    plot(wl,ave_DPF_BB_diff,'LineWidth',3)
    hold on
    plot(wl,DPF_assumed_cyril,'LineWidth',3)
    plot(wl,DPF_homer,'LineWidth',3)
    plot(wl,DPF_scholk_age1,'LineWidth',3)
    plot(wl,DPF_scholk_age2,'LineWidth',3)
    ax = gca;
    ax.FontSize = 20;
    %ylim([0 35]);
    legend('Diffusion value','Broadband toolboxes','Constant value toolboxes','Equation-based, age 20','Equation-based, age 40','Location','Best');
    xlabel('Wavelength (nm)');
    ylabel('Differential Pathlength Factor (DPF)');
    title('Comparison of DPF values used in different toolboxes');
    saveas(gcf,'Comparison of DPF values used in different toolboxes.png')
    
    
    xvector = categorical(["Diffusion value" "Broadband toolboxes" "Constant value toolboxes" "Equation-based, age 20" "Equation-based, age 40"]);
    yvector = [perc_diff_BB_diff(1) perc_diff_BB_cyrilDPF(1) perc_diff_BB_homer(1) perc_diff_BB_scholk1(1) perc_diff_BB_scholk2(1)];
    yvector = abs(yvector);
    
    figure
    bar(xvector,yvector)

%     % Introduce a 5% error in DPF
%     percent_vector = [0.01; 5; 10; 15];
% 
%     for r = 1:length(ave_DPF_BB_diff)
%         five_perc_BB_DPF(r) = (ave_DPF_BB_diff(r)/100) *5;
%         ten_perc_BB_DPF(r) = (ave_DPF_BB_diff(r)/100) *10;
%         fifteen_perc_BB_DPF(r) = (ave_DPF_BB_diff(r)/100) *15;
% 
%         BB_DPF_add5(r) = ave_DPF_BB_diff(r) + five_perc_BB_DPF(r);
%         BB_DPF_add10(r) = ave_DPF_BB_diff(r) + ten_perc_BB_DPF(r);
%         BB_DPF_add15(r) = ave_DPF_BB_diff(r) + fifteen_perc_BB_DPF(r);
% 
%         BB_DPF_take5(r) = ave_DPF_BB_diff(r) - five_perc_BB_DPF(r);
%         BB_DPF_take10(r) = ave_DPF_BB_diff(r) - ten_perc_BB_DPF(r);
%         BB_DPF_take15(r) = ave_DPF_BB_diff(r) - fifteen_perc_BB_DPF(r);
%     end
% 
%     conc_BB_add5percDPF =  E_prime_invert_BB * (-change_I_BB_diff ./ (SD*BB_DPF_add5))';
%     conc_BB_add10percDPF =  E_prime_invert_BB * (-change_I_BB_diff ./ (SD*BB_DPF_add10))';
%     conc_BB_add15percDPF =  E_prime_invert_BB * (-change_I_BB_diff ./ (SD*BB_DPF_add15))';
% 
%     conc_BB_take5percDPF =  E_prime_invert_BB * (-change_I_BB_diff ./ (SD*BB_DPF_take5))';
%     conc_BB_take10percDPF =  E_prime_invert_BB * (-change_I_BB_diff ./ (SD*BB_DPF_take10))';
%     conc_BB_take15percDPF =  E_prime_invert_BB * (-change_I_BB_diff ./ (SD*BB_DPF_take15))';
% 
%     for t = 1:length(ac_conc_change)
%         perc_error_BB_add5DPF(t) = ((conc_BB_add5percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_add10DPF(t) = ((conc_BB_add10percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_add15DPF(t) = ((conc_BB_add15percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_take5DPF(t) = ((conc_BB_take5percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_take10DPF(t) = ((conc_BB_take10percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_take15DPF(t) = ((conc_BB_take15percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPF5(t) = max(perc_error_BB_add5DPF(t),perc_error_BB_take5DPF(t));
%         perc_error_BB_DPF10(t) = max(perc_error_BB_add10DPF(t),perc_error_BB_take10DPF(t));
%         perc_error_BB_DPF15(t) = max(perc_error_BB_add15DPF(t),perc_error_BB_take15DPF(t));
%     end
% 
%     perc_error_DPF_vector = [perc_diff_BB; perc_error_BB_DPF5; perc_error_BB_DPF10; perc_error_BB_DPF15];
%     perc_error_DPF_vector = abs(perc_error_DPF_vector);
% 
%     % perc_error_DPF_vector{l} = perc_error_DPF_vector; % Save each iteration to a cell
%     if l == 1
%         FigH = figure('Position', get(0, 'Screensize'));
%         F = getframe(FigH);
%         imwrite(F.cdata, 'Comparison of percentage errors for increasing DPF percentage changes for Broadband System 1.png', 'png')
%         plot(percent_vector,perc_error_DPF_vector(:,1),'r*-','MarkerSize',20)
%         hold on
%         plot(percent_vector,perc_error_DPF_vector(:,2),'bs-','MarkerSize',10)
%         plot(percent_vector,perc_error_DPF_vector(:,3),'go-','MarkerSize',10)
%         ax = gca;
%         ax.FontSize = 20;
%         ylim([0 35]);
%         legend('HbO','HbR','oxCCO','Location','Best');
%         xlabel('Percentage change in DPF value');
%         ylabel('Percentage error in chromophore concentration values');
%         title({'Comparison of percentage errors for increasing DPF percentage changes', 'for Broadband System 1 (680-921nm)'});
%         saveas(gcf,'Comparison of percentage errors for increasing DPF percentage changes for Broadband System 1.png')
%  
%         for r = 1:length(E)
%             for q = 1:length(ac_conc_change)
%                 five_perc_BB_E(r,q) = (E(r,q)/100) *5;
%                 ten_perc_BB_E(r,q) = (E(r,q)/100) *10;
%                 fifteen_perc_BB_E(r,q) = (E(r,q)/100) *15;
% 
%                 BB_E_add5(r,q) = E(r,q) + five_perc_BB_E(r,q);
%                 BB_E_add10(r,q) = E(r,q) + ten_perc_BB_E(r,q);
%                 BB_E_add15(r,q) = E(r,q) + fifteen_perc_BB_E(r,q);
% 
%                 BB_E_take5(r,q) = E(r,q) - five_perc_BB_E(r,q);
%                 BB_E_take10(r,q) = E(r,q) - ten_perc_BB_E(r,q);
%                 BB_E_take15(r,q) = E(r,q) - fifteen_perc_BB_E(r,q);
%             end
%         end
% 
%             conc_BB_add5percE =  pinv(BB_E_add5) * (-change_I_BB_diff ./ (SD*ave_DPF_BB_diff))';
%             conc_BB_add10percE =  pinv(BB_E_add10) * (-change_I_BB_diff ./ (SD*ave_DPF_BB_diff))';
%             conc_BB_add15percE =  pinv(BB_E_add15) * (-change_I_BB_diff ./ (SD*ave_DPF_BB_diff))';
% 
%             conc_BB_take5percE =  pinv(BB_E_take5) * (-change_I_BB_diff ./ (SD*ave_DPF_BB_diff))';
%             conc_BB_take10percE =  pinv(BB_E_take10) * (-change_I_BB_diff ./ (SD*ave_DPF_BB_diff))';
%             conc_BB_take15percE =  pinv(BB_E_take15) * (-change_I_BB_diff ./ (SD*ave_DPF_BB_diff))';
% 
%         for t = 1:length(ac_conc_change)
%             perc_error_BB_add5E(t) = ((conc_BB_add5percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%             perc_error_BB_add10E(t) = ((conc_BB_add10percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%             perc_error_BB_add15E(t) = ((conc_BB_add15percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%             perc_error_BB_take5E(t) = ((conc_BB_take5percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%             perc_error_BB_take10E(t) = ((conc_BB_take10percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%             perc_error_BB_take15E(t) = ((conc_BB_take15percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%             perc_error_BB_E5(t) = max(perc_error_BB_add5E(t),perc_error_BB_take5E(t));
%             perc_error_BB_E10(t) = max(perc_error_BB_add10E(t),perc_error_BB_take10E(t));
%             perc_error_BB_E15(t) = max(perc_error_BB_add15E(t),perc_error_BB_take15E(t));
%         end
% 
%         perc_error_E_vector_BB = [perc_diff_BB; perc_error_BB_E5; perc_error_BB_E10; perc_error_BB_E15];
%         perc_error_E_vector_BB = abs(perc_error_E_vector_BB);
% 
%         FigH = figure('Position', get(0, 'Screensize'));
%         F = getframe(FigH);
%         imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1.png', 'png')
%         plot(percent_vector,perc_error_E_vector_BB(:,1),'r*-','MarkerSize',20)
%         hold on
%         plot(percent_vector,perc_error_E_vector_BB(:,2),'bs-','MarkerSize',10)
%         plot(percent_vector,perc_error_E_vector_BB(:,3),'go-','MarkerSize',10)
%         ax = gca;
%         ax.FontSize = 20;
%         ylim([0 35]);
%         legend('HbO','HbR','oxCCO','Location','Best');
%         xlabel('Percentage change in extinction coefficient value');
%         ylabel('Percentage error in chromophore concentration values');
%         title({'Comparison of percentage errors for increasing extinction coefficient', 'percentage changes for broadband system 1 (680-921nm)'});
%         saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1 SNR.png')
%     end
% %% Repeat with previously assumed DPF values instead


    %%
% 
%     conc_assDPF_Eadd10 = pinv(BB_E_add10) * (-change_I_BB_diff ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Eadd15_SNR30 = pinv(BB_E_add15) * (-change_I_BB_SNR30 ./ (SD*DPF_assumed))';
% 
%     %
%     conc_DPFadd5_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
% 
%     conc_DPFadd10_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
% 
%     conc_DPFadd15_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
% 
%     %
%     conc_DPFadd10_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
% 
%     conc_DPFadd15_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
% 
%     %
%     %------------------------------------------------------------
%     conc_assDPF_Etake10 = pinv(BB_E_take10) * (-change_I_BB_diff ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
% 
%     conc_assDPF_Etake15_SNR30 = pinv(BB_E_take15) * (-change_I_BB_SNR30 ./ (SD*DPF_assumed))';
% 
%     %
%     conc_DPFtake5_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
% 
%     conc_DPFtake10_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
% 
%     conc_DPFtake15_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
% 
%     %
%     conc_DPFtake10_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
% 
%     conc_DPFtake15_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
% 
%     %
%     %--------------------------------------------------------------------
%     %
%     conc_DPFadd5_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
% 
%     conc_DPFadd10_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
% 
%     conc_DPFadd15_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
% 
%     %
%     conc_DPFadd10_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
% 
%     conc_DPFadd15_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
% 
%     %-------------------------------------------------------------
%     %
%     conc_DPFtake5_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
% 
%     conc_DPFtake10_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
% 
%     conc_DPFtake15_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
% 
%     %
%     conc_DPFtake10_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
% 
%     conc_DPFtake15_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
% 
%     conc_DPFtake15_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
% 
%     conc_DPFtake10_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
% 
%     conc_DPFadd5_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
% 
%     conc_DPFadd5_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
% 
%     conc_DPFtake5_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
% 
%     conc_DPFtake5_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
% 
%     conc_DPFadd10_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
% 
%     conc_DPFadd15_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
% 
%     for t = 1:length(ac_conc_change)
%         perc_error_BB_DPFtake15_Eadd15_SNR40(t) = ((conc_DPFtake15_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFtake10_Eadd10_SNR40(t) = ((conc_DPFtake10_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFtake5_Eadd5_SNR40(t) = ((conc_DPFtake5_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFadd15_Etake15_SNR40(t) = ((conc_DPFadd15_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFadd10_Etake10_SNR40(t) = ((conc_DPFadd10_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFadd5_Etake5_SNR40(t) = ((conc_DPFadd5_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFadd15_Eadd15_SNR40(t) = ((conc_DPFadd15_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFadd10_Eadd10_SNR40(t) = ((conc_DPFadd10_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFadd5_Eadd5_SNR40(t) = ((conc_DPFadd5_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFtake15_Etake15_SNR40(t) = ((conc_DPFtake15_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFtake10_Etake10_SNR40(t) = ((conc_DPFtake10_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFtake5_Etake5_SNR40(t) = ((conc_DPFtake5_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFtake15_Eadd10_SNR40(t) = ((conc_DPFtake15_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFtake15_Eadd5_SNR40(t) = ((conc_DPFtake15_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFtake10_Eadd15_SNR40(t) = ((conc_DPFtake10_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFtake10_Eadd5_SNR40(t) = ((conc_DPFtake10_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFadd5_Etake15_SNR40(t) = ((conc_DPFadd5_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFadd5_Etake10_SNR40(t) = ((conc_DPFadd5_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFtake5_Eadd15_SNR40(t) = ((conc_DPFtake5_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFtake5_Eadd10_SNR40(t) = ((conc_DPFtake5_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFadd10_Etake15_SNR40(t) = ((conc_DPFadd10_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFadd10_Etake5_SNR40(t) = ((conc_DPFadd10_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_DPFadd15_Etake10_SNR40(t) = ((conc_DPFadd15_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
%         perc_error_BB_DPFadd15_Etake5_SNR40(t) = ((conc_DPFadd15_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
% 
%         perc_error_BB_5s_SNR40(t) = max(abs([perc_error_BB_DPFtake5_Etake5_SNR40(t), perc_error_BB_DPFadd5_Eadd5_SNR40(t), perc_error_BB_DPFadd5_Etake5_SNR40(t), perc_error_BB_DPFtake5_Eadd5_SNR40(t)]));
%         perc_error_BB_10s_SNR40(t) = max(abs([perc_error_BB_DPFtake10_Etake10_SNR40(t), perc_error_BB_DPFadd10_Eadd10_SNR40(t), perc_error_BB_DPFadd10_Etake10_SNR40(t), perc_error_BB_DPFtake10_Eadd10_SNR40(t)]));
%         perc_error_BB_15s_SNR40(t) = max(abs([perc_error_BB_DPFtake15_Etake15_SNR40(t), perc_error_BB_DPFadd15_Eadd15_SNR40(t), perc_error_BB_DPFadd15_Etake15_SNR40(t), perc_error_BB_DPFtake15_Eadd15_SNR40(t)]));
%     end
%     %%
%     perc_accum_errors_DPFE_BB = [perc_diff_BB; perc_error_BB_5s_SNR40; perc_error_BB_10s_SNR40; perc_error_BB_15s_SNR40];
% 
%     perc_accum_errors_DPFE_BB_saved{l} = perc_accum_errors_DPFE_BB(:,:); % Save iteration to a cell array
% 
%     % FigH = figure('Position', get(0, 'Screensize'));
%     % F = getframe(FigH);
%     % imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient and DPF percentage changes for broadband system 1 SNR 40.png', 'png')
%     % plot(percent_vector,perc_accum_errors_DPFE_BB(:,1),'r*-','MarkerSize',20)
%     % hold on
%     % plot(percent_vector,perc_accum_errors_DPFE_BB(:,2),'bs-','MarkerSize',10)
%     % plot(percent_vector,perc_accum_errors_DPFE_BB(:,5),'go-','MarkerSize',10)
%     % ax = gca;
%     % ax.FontSize = 20;
%     % %ylim([0 35]);
%     % legend('HbO','HbR','oxCCO','Location','Best');
%     % xlabel('Percentage change in both extinction coefficient and DPF values');
%     % ylabel('Percentage error in chromophore concentration values');
%     % title({'Comparison of percentage errors for increasing extinction coefficient', 'and DPF percentage changes for broadband system 1 (680-921nm), SNR = 40'});
%     % saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1 SNR 40.png')

    %% Again for SNR 50


    conc_assDPF_Eadd10 = pinv(BB_E_add10) * (-change_I_BB_diff ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Eadd5_SNR50 = pinv(BB_E_add5) * (-change_I_BB_SNR50 ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Eadd10_SNR50 = pinv(BB_E_add10) * (-change_I_BB_SNR50 ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Eadd15_SNR50 = pinv(BB_E_add15) * (-change_I_BB_SNR50 ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Eadd15_SNR30 = pinv(BB_E_add15) * (-change_I_BB_SNR30 ./ (SD*DPF_assumed_cyril))';

    %
    conc_DPFadd5_Eadd5_SNR50 = pinv(BB_E_add5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add5))';

    conc_DPFadd10_Eadd5_SNR50 = pinv(BB_E_add5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add10))';

    conc_DPFadd15_Eadd5_SNR50 = pinv(BB_E_add5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add15))';

    %
    conc_DPFadd10_Eadd10_SNR50 = pinv(BB_E_add10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add10))';

    conc_DPFadd15_Eadd15_SNR50 = pinv(BB_E_add15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add15))';

    %
    %------------------------------------------------------------
    conc_assDPF_Etake10 = pinv(BB_E_take10) * (-change_I_BB_diff ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Etake5_SNR50 = pinv(BB_E_take5) * (-change_I_BB_SNR50 ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Etake10_SNR50 = pinv(BB_E_take10) * (-change_I_BB_SNR50 ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Etake15_SNR50 = pinv(BB_E_take15) * (-change_I_BB_SNR50 ./ (SD*DPF_assumed_cyril))';

    conc_assDPF_Etake15_SNR30 = pinv(BB_E_take15) * (-change_I_BB_SNR30 ./ (SD*DPF_assumed_cyril))';

    %
    conc_DPFtake5_Etake5_SNR50 = pinv(BB_E_take5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take5))';

    conc_DPFtake10_Etake5_SNR50 = pinv(BB_E_take5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take10))';

    conc_DPFtake15_Etake5_SNR50 = pinv(BB_E_take5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take15))';

    %
    conc_DPFtake10_Etake10_SNR50 = pinv(BB_E_take10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take10))';

    conc_DPFtake15_Etake15_SNR50 = pinv(BB_E_take15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take15))';

    %
    %--------------------------------------------------------------------
    %
    conc_DPFadd5_Etake5_SNR50 = pinv(BB_E_take5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add5))';

    conc_DPFadd10_Etake5_SNR50 = pinv(BB_E_take5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add10))';

    conc_DPFadd15_Etake5_SNR50 = pinv(BB_E_take5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add15))';

    %
    conc_DPFadd10_Etake10_SNR50 = pinv(BB_E_take10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add10))';

    conc_DPFadd15_Etake15_SNR50 = pinv(BB_E_take15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add15))';

    %-------------------------------------------------------------
    %
    conc_DPFtake5_Eadd5_SNR50 = pinv(BB_E_add5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take5))';

    conc_DPFtake10_Eadd5_SNR50 = pinv(BB_E_add5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take10))';

    conc_DPFtake15_Eadd5_SNR50 = pinv(BB_E_add5) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take15))';

    %
    conc_DPFtake10_Eadd10_SNR50 = pinv(BB_E_add10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take10))';

    conc_DPFtake15_Eadd15_SNR50 = pinv(BB_E_add15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take15))';

    conc_DPFtake15_Eadd10_SNR50 = pinv(BB_E_add10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take15))';

    conc_DPFtake10_Eadd15_SNR50 = pinv(BB_E_add15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take10))';

    conc_DPFadd5_Etake15_SNR50 = pinv(BB_E_take15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add5))';

    conc_DPFadd5_Etake10_SNR50 = pinv(BB_E_take10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add5))';

    conc_DPFtake5_Eadd15_SNR50 = pinv(BB_E_add15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take5))';

    conc_DPFtake5_Eadd10_SNR50 = pinv(BB_E_add10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_take5))';

    conc_DPFadd10_Etake15_SNR50 = pinv(BB_E_take15) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add10))';

    conc_DPFadd15_Etake10_SNR50 = pinv(BB_E_take10) * (-change_I_BB_SNR50 ./ (SD*BB_DPF_add15))';

    for t = 1:length(ac_conc_change)
        perc_error_BB_DPFtake15_Eadd15_SNR50(t) = ((conc_DPFtake15_Eadd15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFtake10_Eadd10_SNR50(t) = ((conc_DPFtake10_Eadd10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFtake5_Eadd5_SNR50(t) = ((conc_DPFtake5_Eadd5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFadd15_Etake15_SNR50(t) = ((conc_DPFadd15_Etake15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFadd10_Etake10_SNR50(t) = ((conc_DPFadd10_Etake10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFadd5_Etake5_SNR50(t) = ((conc_DPFadd5_Etake5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFadd15_Eadd15_SNR50(t) = ((conc_DPFadd15_Eadd15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFadd10_Eadd10_SNR50(t) = ((conc_DPFadd10_Eadd10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFadd5_Eadd5_SNR50(t) = ((conc_DPFadd5_Eadd5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFtake15_Etake15_SNR50(t) = ((conc_DPFtake15_Etake15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFtake10_Etake10_SNR50(t) = ((conc_DPFtake10_Etake10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFtake5_Etake5_SNR50(t) = ((conc_DPFtake5_Etake5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFtake15_Eadd10_SNR50(t) = ((conc_DPFtake15_Eadd10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFtake15_Eadd5_SNR50(t) = ((conc_DPFtake15_Eadd5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFtake10_Eadd15_SNR50(t) = ((conc_DPFtake10_Eadd15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFtake10_Eadd5_SNR50(t) = ((conc_DPFtake10_Eadd5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFadd5_Etake15_SNR50(t) = ((conc_DPFadd5_Etake15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFadd5_Etake10_SNR50(t) = ((conc_DPFadd5_Etake10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFtake5_Eadd15_SNR50(t) = ((conc_DPFtake5_Eadd15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFtake5_Eadd10_SNR50(t) = ((conc_DPFtake5_Eadd10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFadd10_Etake15_SNR50(t) = ((conc_DPFadd10_Etake15_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFadd10_Etake5_SNR50(t) = ((conc_DPFadd10_Etake5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_DPFadd15_Etake10_SNR50(t) = ((conc_DPFadd15_Etake10_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        perc_error_BB_DPFadd15_Etake5_SNR50(t) = ((conc_DPFadd15_Etake5_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;

        perc_error_BB_5s_SNR50(t) = max(abs([perc_error_BB_DPFtake5_Etake5_SNR50(t), perc_error_BB_DPFadd5_Eadd5_SNR50(t), perc_error_BB_DPFadd5_Etake5_SNR50(t), perc_error_BB_DPFtake5_Eadd5_SNR50(t)]));
        perc_error_BB_10s_SNR50(t) = max(abs([perc_error_BB_DPFtake10_Etake10_SNR50(t), perc_error_BB_DPFadd10_Eadd10_SNR50(t), perc_error_BB_DPFadd10_Etake10_SNR50(t), perc_error_BB_DPFtake10_Eadd10_SNR50(t)]));
        perc_error_BB_15s_SNR50(t) = max(abs([perc_error_BB_DPFtake15_Etake15_SNR50(t), perc_error_BB_DPFadd15_Eadd15_SNR50(t), perc_error_BB_DPFadd15_Etake15_SNR50(t), perc_error_BB_DPFtake15_Eadd15_SNR50(t)]));
    end
    %%
    perc_accum_errors_DPFE_BBSNR50 = [perc_diff_BB_diff; perc_error_BB_5s_SNR50; perc_error_BB_10s_SNR50; perc_error_BB_15s_SNR50];

    perc_accum_errors_DPFE_BBSNR50_saved{l} = perc_accum_errors_DPFE_BBSNR50(:,:); %Save each iteration to cell array

    % FigH = figure('Position', get(0, 'Screensize'));
    % F = getframe(FigH);
    % imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient and DPF percentage changes for broadband system 1 SNR 50.png', 'png')
    % plot(percent_vector,perc_accum_errors_DPFE_BBSNR50(:,1),'r*-','MarkerSize',20)
    % hold on
    % plot(percent_vector,perc_accum_errors_DPFE_BBSNR50(:,2),'bs-','MarkerSize',10)
    % plot(percent_vector,perc_accum_errors_DPFE_BBSNR50(:,5),'go-','MarkerSize',10)
    % ax = gca;
    % ax.FontSize = 20;
    % %ylim([0 35]);
    % legend('HbO','HbR','oxCCO','Location','Best');
    % xlabel('Percentage change in both extinction coefficient and DPF values');
    % ylabel('Percentage error in chromophore concentration values');
    % title({'Comparison of percentage errors for increasing extinction coefficient', 'and DPF percentage changes for broadband system 1 (680-921nm), SNR = 50'});
    % saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1 SNR 50.png')

end

%%

for u = 1:4
    for j = 1:3
        for k = 1:10
            vector_error_BBSNR50(k,j) = perc_accum_errors_DPFE_BBSNR50_saved{1,k}(u,j);
        end
        ave_vector_error_BBSNR50(u,j) = mean(vector_error_BBSNR50(:,j));
        high_bar_vector_error_BBSNR50(u,j) = max(vector_error_BBSNR50(:,j));
        low_bar_vector_error_BBSNR50(u,j) = min(vector_error_BBSNR50(:,j));
    end
end

%%
FigH = figure('Position', get(0, 'Screensize'));
F = getframe(FigH);
imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient and DPF percentage changes for broadband system 1 SNR 50.png', 'png')
errorbar(percent_vector,ave_vector_error_BBSNR50(:,1),low_bar_vector_error_BBSNR50(:,1),high_bar_vector_error_BBSNR50(:,1),'r*-','MarkerSize',20,'LineWidth',1.8)
hold on
errorbar(percent_vector,ave_vector_error_BBSNR50(:,2),low_bar_vector_error_BBSNR50(:,2),high_bar_vector_error_BBSNR50(:,2),'bs-','MarkerSize',20,'LineWidth',1.5)
errorbar(percent_vector,ave_vector_error_BBSNR50(:,3),low_bar_vector_error_BBSNR50(:,3),high_bar_vector_error_BBSNR50(:,3),'go-','MarkerSize',20,'LineWidth',1)
grid on
ax = gca;
ax.FontSize = 20;
xlim([-5 20]);
legend('HbO','HbR','oxCCO','Location','Best');
xlabel('Percentage change in both extinction coefficient and DPF values');
ylabel('Percentage error in chromophore concentration values');
title({'Comparison of percentage errors for increasing extinction coefficient', 'and DPF percentage changes for broadband system 1 (680-921nm), SNR = 50'});
saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1 SNR 50.png')
