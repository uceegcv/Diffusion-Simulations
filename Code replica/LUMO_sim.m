%% Diffuse reflection model for a slab geometry
% % September 2023
% clear all;
% clc;
% close all;
% addpath(genpath('N:\Matlab code\homer2_v2_8_11022018')) 


%% Code for Broadband System 1

for nwavs = 2:200
    for pert_number = 1:100

        % Enter simulation data
        SD = 30; % Source-detector separation on same surface of slab (mm)
        slab = 50; % Slab thickness (mm)
        % x contains all broadband wavelengths
        x = 680:921;
        
        % Contains broadband wavelengths of interest - the one to edit
        %wl = 680:921; %nm (cannot be a larger range than x)
        wl = 680:921; % de Roever
        
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
        
        E_rewrite = [E(:,1) E(:,2)]; % ADD
        E = E_rewrite;
        
        E_new_rewrite = [E_new(:,1) E_new(:,2)]; % ADD
        E_new = E_new_rewrite;
        
        % Extinction coefficients
        ext_coeffs = tissue_specific_extinction_coefficient_650to1000;
        
        % water, HbO2, HHb, CCO
        ext_coeffs_CYRIL = ext_coeffs(wl(1)-649:wl(end)-649,2:5);
        ext_coeffs_CYRIL_fat = ext_coeffs(wl(1)-649:wl(end)-649,8);
        
        % Re-order
        ext_CYRIL = [ext_coeffs_CYRIL(:,2) ext_coeffs_CYRIL(:,3) ext_coeffs_CYRIL(:,1) ext_coeffs_CYRIL_fat ext_coeffs_CYRIL(:,4)];
        ext_CYRIL(:,1:2) = ext_CYRIL(:,1:2)/(10e5*2*2.3);
        ext_CYRIL(:,4:5) = ext_CYRIL(:,4:5)/(10e5*2*2.3);
        ext_CYRIL(:,3) = ext_CYRIL(:,3);
        
        % ext_BioBook = tissue_specific_extinction_coefficient_650to1000_OpticsBook;
        % %%
        % for i = 1:5
        %     %ext_BioBook_extend = interp1(ext_BioBook(:,1),ext_BioBook(:,i+1),ext_BioBook(1,1):ext_BioBook(end,1));
        %     ext_BioBook_extend = interp1(ext_BioBook(:,1),ext_BioBook(:,i+1),680:921);
        % 
        %     ext_BioBook_new(:,i) = ext_BioBook_extend.';
        % end
        % 
        % ext_BioBook_new(:,1:2) = ext_BioBook_new(:,1:2)/10e6;
        % ext_BioBook_new(:,3) = ext_BioBook_new(:,3)/1000;
        % ext_BioBook_new(:,4) = ext_BioBook_new(:,4)/100;
        % ext_BioBook_new(:,5) = ext_BioBook_new(:,5)/10e6;
        
        %ext_CYRIL = ext_CYRIL/(2*2.3);
        % %% Chromophore spectra
        %Plot constituent spectra
    %     FigH = figure('Position', get(0, 'Screensize'));
    %     F = getframe(FigH);
    %     imwrite(F.cdata, 'Constituent spectra (and polynomial fits) from Homer2 for broadband system 1.png', 'png')
    %     plot(x,E(:,1),'ro-');hold on;
    %     plot(x,E(:,2),'bo-');
    %     % plot(x,(E(:,3)/10),'go-'); %Divide by 10 to improve display
    %     % plot(x,(E(:,4)/10),'ko-'); %Divide by 10 to improve display
    %     % plot(x,E(:,3),'mo-');
    %     xlim([min(x) max(x)])
    %     legend('HbO','HbR');
    %     xlabel('Wavelength (nm)');
    %     ylabel('Specific abs. coeff or Abs. coeff/10 (mm^{-1}M^{-1} or mm^{-1})');
    %     title('Constituent spectra from Homer2 for broadband system 1 (680-921nm)');
    %     saveas(gcf,'Constituent spectra (and polynomial fits) from Homer2 for broadband system 1.png')
    %     
        % 
        % figure;
        % plot(x,ext_CYRIL(:,1),'ro-');hold on;
        % plot(x,ext_CYRIL(:,2),'bo-');
        % plot(x,(ext_CYRIL(:,3)),'go-'); %Divide by 10 to improve display
        % plot(x,(ext_CYRIL(:,4)),'ko-'); %Divide by 10 to improve display
        % plot(x,ext_CYRIL(:,5),'mo-');
        % xlim([min(x) max(x)])
        % legend('HbO','HbR','H2O','Lipid','aa3');
        % xlabel('Wavelength (nm)');
        % ylabel('Specific abs. coeff or Abs. coeff/10 (mm^{-1}M^{-1} or mm^{-1})');
        % title('Constituent spectra (and polynomial fits) from CYRIL function');
        % 
        % figure;
        % plot(x,ext_CYRIL(:,1),'ro-');hold on;
        % plot(x,ext_CYRIL(:,2),'bo-');
        % plot(x,ext_CYRIL(:,5),'mo-');
        % plot(x,E(:,1),'ro-');
        % plot(x,E(:,2),'bo-');
        % plot(x,E(:,5),'mo-');
        % xlim([min(x) max(x)])
        % legend('HbO','HbR','aa3','HbO','HbR','aa3');
        % xlabel('Wavelength (nm)');
        % ylabel('Specific abs. coeff or Abs. coeff/10 (mm^{-1}M^{-1} or mm^{-1})');
        % title('Constituent spectra comparing Homer2 and CYRIL values');
        % 
        % figure;
        % plot(x,ext_CYRIL(:,1),'ro-');hold on;
        % plot(x,ext_CYRIL(:,2),'bo-');
        % plot(x,ext_CYRIL(:,5),'mo-');
        % plot(x,E(:,1),'rs-');
        % plot(x,E(:,2),'bs-');
        % plot(x,E(:,5),'ms-');
        % plot(x,ext_BioBook_new(:,1),'r*-');
        % plot(x,ext_BioBook_new(:,2),'b*-');
        % plot(x,ext_BioBook_new(:,5),'m*-');
        % xlim([min(x) max(x)])
        % legend('HbO CYRIL','HbR CYRIL','aa3 CYRIL','HbO Homer2','HbR Homer2','aa3 Homer2','HbO Book','HbR Book','aa3 Book');
        % xlabel('Wavelength (nm)');
        % ylabel('Specific abs. coeff or Abs. coeff/10 (mm^{-1}M^{-1} or mm^{-1})');
        % title('Constituent spectra comparing Homer2, CYRIL and Book values');
        
        %% Code for Broadband System
        
        W(1) = 0.8;    %Water fraction
        L1(1) = 0.116;  %Lipid fraction
        B(1) = 0.012;  %Background non-wavelength-dependent absorption coefficient (mm-1)
        C_HbO(1) = 56; %Concentration of HbO (uM)
        C_HbR(1) = 24; %Concentration of HbR (uM)
        C_aa3(1) = 4.9; %Concentration of aa3 (uM)
        % mua1 = (E(:,1)*C_HbO(1)) + (E(:,2)*C_HbR(1)) + (E(:,3)*W(1)) + (E(:,4)*L1(1)) + (E(:,5)*C_aa3(1));
        mua1 = (E(:,1)*C_HbO(1)) + (E(:,2)*C_HbR(1));
        
        % Put concs in vector
        % concs1 = [C_HbO(1) C_HbR(1) W(1) L1(1) C_aa3(1)];
        concs1 = [C_HbO(1) C_HbR(1)];
        
        I_original = exp(-mua1);
        
        %% Change concentrations
            
       % Induce random changes
        % Change HbO between +1.5 and -1.5 uM
        HbO_pert = 0.01:0.01:3.01; % define the numbers
        HbO_change = 1.51 - HbO_pert(randi([1,numel(HbO_pert)]));
        
        % Change HbR between +1 and -1 uM
        HbR_pert = 0.01:0.01:2.01; % define the numbers
        HbR_change = 1.01 - HbR_pert(randi([1,numel(HbR_pert)]));
        
        % Change CCO between +0.3 and -0.3 uM
        CCO_pert = 0.01:0.01:0.61; % define the numbers
        CCO_change = 0.31 - CCO_pert(randi([1,numel(CCO_pert)]));
    
        W(2) = 0.8;    %Water fraction
        L1(2) = 0.116;  %Lipid fraction
        B(2) = 0.012;  %Background non-wavelength-dependent absorption coefficient (mm-1)
        C_HbO(2) = C_HbO(1) + HbO_change; %Concentration of HbO (uM)
        C_HbR(2) = C_HbR(1) + HbR_change; %Concentration of HbR (uM)
        C_aa3(2) = C_aa3(1) + CCO_change; %Concentration of aa3 (uM)
        % mua2 = (E(:,1)*C_HbO(2)) + (E(:,2)*C_HbR(2)) + (E(:,3)*W(2)) + (E(:,4)*L1(2)) + (E(:,5)*C_aa3(2));
        mua2 = (E(:,1)*C_HbO(2)) + (E(:,2)*C_HbR(2));
        
        % concs2 = [C_HbO(2) C_HbR(2) W(2) L1(2) C_aa3(2)];
        concs2 = [C_HbO(2) C_HbR(2)];
        
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
        
        change_I_BB = log(R2_BB) - log(R1_BB);
        
        ave_DPF_BB = (DPF1_BB + DPF2_BB)/2;
        
        for t = 1:length(x)
            calc_I_BB(t) = ave_DPF_BB(t)*SD*(C_HbO(2)-C_HbO(1))*E(t,1);
        end
        
        for t = 1:length(x)
            perc_diff_I(t) = 100 - ((calc_I_BB(t)/change_I_BB(t))*100);
        end
        
        E_prime_invert_BB = pinv(E);
        
        % This is the value to pay attention to!! From diffusion
        conc_BB =  E_prime_invert_BB * (-change_I_BB ./ (SD*ave_DPF_BB))'; % From diffusion
        
        %%
        for t = 1:length(ac_conc_change)
            perc_diff_BB(t) = ((conc_BB(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        end
        
        
        SNR_vector = [20; 30; 40; 50; 60];
        
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
        
        ave_DPF_BB = (DPF1_BB + DPF2_BB)/2;
        
        for t = 1:length(x)
            calc_I_BB(t) = ave_DPF_BB(t)*SD*(C_HbO(2)-C_HbO(1))*E(t,1);
        end
        
        for t = 1:length(x)
            perc_diff_I_BB(t) = 100 - ((calc_I_BB(t)/change_I_BB(t))*100);
        end
        
        % This is the value to pay attention to!! From diffusion
        conc_BB =  E_prime_invert_BB * (-change_I_BB ./ (SD*ave_DPF_BB))'; % From diffusion
        
        conc_BB_SNR20 =  E_prime_invert_BB * (-change_I_BB_SNR20 ./ (SD*ave_DPF_BB))';
        conc_BB_SNR30 =  E_prime_invert_BB * (-change_I_BB_SNR30 ./ (SD*ave_DPF_BB))';
        conc_BB_SNR40 =  E_prime_invert_BB * (-change_I_BB_SNR40 ./ (SD*ave_DPF_BB))';
        conc_BB_SNR50 =  E_prime_invert_BB * (-change_I_BB_SNR50 ./ (SD*ave_DPF_BB))';
        conc_BB_SNR60 =  E_prime_invert_BB * (-change_I_BB_SNR60 ./ (SD*ave_DPF_BB))';
        
        
        for t = 1:length(ac_conc_change)
            perc_diff_BB(t) = ((conc_BB(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_diff_BB_SNR20(t) = ((conc_BB_SNR20(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_diff_BB_SNR30(t) = ((conc_BB_SNR30(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_diff_BB_SNR40(t) = ((conc_BB_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_diff_BB_SNR50(t) = ((conc_BB_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_diff_BB_SNR60(t) = ((conc_BB_SNR60(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        end
        
        perc_vector_SNR = [20;30;40;50;60;100];
        perc_diff_SNR_vector = [perc_diff_BB_SNR20;perc_diff_BB_SNR30;perc_diff_BB_SNR40;perc_diff_BB_SNR50;perc_diff_BB_SNR60;perc_diff_BB];
        perc_diff_SNR_vector = abs(perc_diff_SNR_vector);
        
    %     FigH = figure('Position', get(0, 'Screensize'));
    %     F = getframe(FigH);
    %     imwrite(F.cdata, 'Comparison of percentage errors for increasing SNR of broadband system 1.png', 'png')
    %     plot(perc_vector_SNR,perc_diff_SNR_vector(:,1),'r*-','MarkerSize',20)
    %     hold on
    %     plot(perc_vector_SNR,perc_diff_SNR_vector(:,2),'bs-','MarkerSize',10)
    %     % plot(perc_vector_SNR,perc_diff_SNR_vector(:,3),'go-','MarkerSize',10)
    %     ax = gca;
    %     ax.FontSize = 20;
    %     legend('HbO','HbR','Location','Best');
    %     xlabel('Absolute SNR value');
    %     ylabel('Percentage error in chromophore concentration values');
    %     title({'Comparison of percentage errors for increasing SNR', 'of broadband system 1 (680-921nm)'});
    %     saveas(gcf,'Comparison of percentage errors for increasing SNR of broadband system 1.png')
        
        %% Introduce errors in BB measurements due to DPF changes
        
        % Introduce a 5% error in DPF
        percent_vector = [0.01; 5; 10; 15];
        
        for r = 1:length(ave_DPF_BB)
            five_perc_BB_DPF(r) = (ave_DPF_BB(r)/100) *5;
            ten_perc_BB_DPF(r) = (ave_DPF_BB(r)/100) *10;
            fifteen_perc_BB_DPF(r) = (ave_DPF_BB(r)/100) *15;
            
            BB_DPF_add5(r) = ave_DPF_BB(r) + five_perc_BB_DPF(r);
            BB_DPF_add10(r) = ave_DPF_BB(r) + ten_perc_BB_DPF(r);
            BB_DPF_add15(r) = ave_DPF_BB(r) + fifteen_perc_BB_DPF(r);
            
            BB_DPF_take5(r) = ave_DPF_BB(r) - five_perc_BB_DPF(r);
            BB_DPF_take10(r) = ave_DPF_BB(r) - ten_perc_BB_DPF(r);
            BB_DPF_take15(r) = ave_DPF_BB(r) - fifteen_perc_BB_DPF(r);
        end
        
        conc_BB_add5percDPF =  E_prime_invert_BB * (-change_I_BB ./ (SD*BB_DPF_add5))';
        conc_BB_add10percDPF =  E_prime_invert_BB * (-change_I_BB ./ (SD*BB_DPF_add10))';
        conc_BB_add15percDPF =  E_prime_invert_BB * (-change_I_BB ./ (SD*BB_DPF_add15))';
        
        conc_BB_take5percDPF =  E_prime_invert_BB * (-change_I_BB ./ (SD*BB_DPF_take5))';
        conc_BB_take10percDPF =  E_prime_invert_BB * (-change_I_BB ./ (SD*BB_DPF_take10))';
        conc_BB_take15percDPF =  E_prime_invert_BB * (-change_I_BB ./ (SD*BB_DPF_take15))';
        
        for t = 1:length(ac_conc_change)
            perc_error_BB_add5DPF(t) = ((conc_BB_add5percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_add10DPF(t) = ((conc_BB_add10percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_add15DPF(t) = ((conc_BB_add15percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_take5DPF(t) = ((conc_BB_take5percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_take10DPF(t) = ((conc_BB_take10percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_take15DPF(t) = ((conc_BB_take15percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPF5(t) = max(perc_error_BB_add5DPF(t),perc_error_BB_take5DPF(t));
            perc_error_BB_DPF10(t) = max(perc_error_BB_add10DPF(t),perc_error_BB_take10DPF(t));
            perc_error_BB_DPF15(t) = max(perc_error_BB_add15DPF(t),perc_error_BB_take15DPF(t));
        end
        
        perc_error_DPF_vector = [perc_diff_BB; perc_error_BB_DPF5; perc_error_BB_DPF10; perc_error_BB_DPF15];
        perc_error_DPF_vector = abs(perc_error_DPF_vector);
        
    %     FigH = figure('Position', get(0, 'Screensize'));
    %     F = getframe(FigH);
    %     imwrite(F.cdata, 'Comparison of percentage errors for increasing DPF percentage changes for Broadband System 1.png', 'png')
    %     plot(percent_vector,perc_error_DPF_vector(:,1),'r*-','MarkerSize',20)
    %     hold on
    %     plot(percent_vector,perc_error_DPF_vector(:,2),'bs-','MarkerSize',10)
    %     % plot(percent_vector,perc_error_DPF_vector(:,3),'go-','MarkerSize',10)
    %     ax = gca;
    %     ax.FontSize = 20;
    %     ylim([0 35]);
    %     legend('HbO','HbR','Location','Best');
    %     xlabel('Percentage change in DPF value');
    %     ylabel('Percentage error in chromophore concentration values');
    %     title({'Comparison of percentage errors for increasing DPF percentage changes', 'for Broadband System 1 (680-921nm)'});
    %     saveas(gcf,'Comparison of percentage errors for increasing DPF percentage changes for Broadband System 1.png')
        
        for r = 1:length(E)
            for q = 1:length(ac_conc_change)
                five_perc_BB_E(r,q) = (E(r,q)/100) *5;
                ten_perc_BB_E(r,q) = (E(r,q)/100) *10;
                fifteen_perc_BB_E(r,q) = (E(r,q)/100) *15;
        
                BB_E_add5(r,q) = E(r,q) + five_perc_BB_E(r,q);
                BB_E_add10(r,q) = E(r,q) + ten_perc_BB_E(r,q);
                BB_E_add15(r,q) = E(r,q) + fifteen_perc_BB_E(r,q);
        
                BB_E_take5(r,q) = E(r,q) - five_perc_BB_E(r,q);
                BB_E_take10(r,q) = E(r,q) - ten_perc_BB_E(r,q);
                BB_E_take15(r,q) = E(r,q) - fifteen_perc_BB_E(r,q);
            end
        end
        
        conc_BB_add5percE =  pinv(BB_E_add5) * (-change_I_BB ./ (SD*ave_DPF_BB))';
        conc_BB_add10percE =  pinv(BB_E_add10) * (-change_I_BB ./ (SD*ave_DPF_BB))';
        conc_BB_add15percE =  pinv(BB_E_add15) * (-change_I_BB ./ (SD*ave_DPF_BB))';
        
        conc_BB_take5percE =  pinv(BB_E_take5) * (-change_I_BB ./ (SD*ave_DPF_BB))';
        conc_BB_take10percE =  pinv(BB_E_take10) * (-change_I_BB ./ (SD*ave_DPF_BB))';
        conc_BB_take15percE =  pinv(BB_E_take15) * (-change_I_BB ./ (SD*ave_DPF_BB))';
        
        for t = 1:length(ac_conc_change)
            perc_error_BB_add5E(t) = ((conc_BB_add5percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_add10E(t) = ((conc_BB_add10percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_add15E(t) = ((conc_BB_add15percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_take5E(t) = ((conc_BB_take5percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_take10E(t) = ((conc_BB_take10percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_take15E(t) = ((conc_BB_take15percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_E5(t) = max(perc_error_BB_add5E(t),perc_error_BB_take5E(t));
            perc_error_BB_E10(t) = max(perc_error_BB_add10E(t),perc_error_BB_take10E(t));
            perc_error_BB_E15(t) = max(perc_error_BB_add15E(t),perc_error_BB_take15E(t));
        end
        
        perc_error_E_vector_BB = [perc_diff_BB; perc_error_BB_E5; perc_error_BB_E10; perc_error_BB_E15];
        perc_error_E_vector_BB = abs(perc_error_E_vector_BB);
        
    %     FigH = figure('Position', get(0, 'Screensize'));
    %     F = getframe(FigH);
    %     imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1.png', 'png')
    %     plot(percent_vector,perc_error_E_vector_BB(:,1),'r*-','MarkerSize',20)
    %     hold on
    %     plot(percent_vector,perc_error_E_vector_BB(:,2),'bs-','MarkerSize',10)
    %     % plot(percent_vector,perc_error_E_vector_BB(:,3),'go-','MarkerSize',10)
    %     ax = gca;
    %     ax.FontSize = 20;
    %     ylim([0 35]);
    %     legend('HbO','HbR','Location','Best');
    %     xlabel('Percentage change in extinction coefficient value');
    %     ylabel('Percentage error in chromophore concentration values');
    %     title({'Comparison of percentage errors for increasing extinction coefficient', 'percentage changes for broadband system 1 (680-921nm)'});
    %     saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1.png')
        
        % %% Repeat with previously assumed DPF values instead
        
        DPF_dep = stretch_DPF_Lambda_Dependency_680to915; 
        
        for r = 1:length(wl)
            DPF_index_BB(r) = find(DPF_dep(:,1) == wl(r).');
            DPF_vals_BB(r) = DPF_dep(DPF_index_BB(r),2);
        end
        
        set_DPF = 4.99;
        
        DPF_assumed = DPF_vals_BB * set_DPF;
        
        conc_assumedDPF_BB = E_prime_invert_BB * (-change_I_BB ./ (SD*DPF_assumed))';
        
        %%
        
        conc_assDPF_Eadd10 = pinv(BB_E_add10) * (-change_I_BB ./ (SD*DPF_assumed))';
        
        conc_assDPF_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
        
        conc_assDPF_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
        
        conc_assDPF_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
        
        conc_assDPF_Eadd15_SNR30 = pinv(BB_E_add15) * (-change_I_BB_SNR30 ./ (SD*DPF_assumed))';
        
        %
        conc_DPFadd5_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
        
        conc_DPFadd10_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
        
        conc_DPFadd15_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
        
        %
        conc_DPFadd10_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
        
        conc_DPFadd15_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
        
        %
        %------------------------------------------------------------
        conc_assDPF_Etake10 = pinv(BB_E_take10) * (-change_I_BB ./ (SD*DPF_assumed))';
        
        conc_assDPF_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
        
        conc_assDPF_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
        
        conc_assDPF_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*DPF_assumed))';
        
        conc_assDPF_Etake15_SNR30 = pinv(BB_E_take15) * (-change_I_BB_SNR30 ./ (SD*DPF_assumed))';
        
        %
        conc_DPFtake5_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
        
        conc_DPFtake10_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
        
        conc_DPFtake15_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
        
        %
        conc_DPFtake10_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
        
        conc_DPFtake15_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
        
        %
        %--------------------------------------------------------------------
        %
        conc_DPFadd5_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
        
        conc_DPFadd10_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
        
        conc_DPFadd15_Etake5_SNR40 = pinv(BB_E_take5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
        
        %
        conc_DPFadd10_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
        
        conc_DPFadd15_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
        
        %-------------------------------------------------------------
        %
        conc_DPFtake5_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
        
        conc_DPFtake10_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
        
        conc_DPFtake15_Eadd5_SNR40 = pinv(BB_E_add5) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
        
        %
        conc_DPFtake10_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
        
        conc_DPFtake15_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
        
        conc_DPFtake15_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take15))';
        
        conc_DPFtake10_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take10))';
        
        conc_DPFadd5_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
        
        conc_DPFadd5_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add5))';
        
        conc_DPFtake5_Eadd15_SNR40 = pinv(BB_E_add15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
        
        conc_DPFtake5_Eadd10_SNR40 = pinv(BB_E_add10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_take5))';
        
        conc_DPFadd10_Etake15_SNR40 = pinv(BB_E_take15) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add10))';
        
        conc_DPFadd15_Etake10_SNR40 = pinv(BB_E_take10) * (-change_I_BB_SNR40 ./ (SD*BB_DPF_add15))';
        
        for t = 1:length(ac_conc_change)
            perc_error_BB_DPFtake15_Eadd15_SNR40(t) = ((conc_DPFtake15_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFtake10_Eadd10_SNR40(t) = ((conc_DPFtake10_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFtake5_Eadd5_SNR40(t) = ((conc_DPFtake5_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFadd15_Etake15_SNR40(t) = ((conc_DPFadd15_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFadd10_Etake10_SNR40(t) = ((conc_DPFadd10_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFadd5_Etake5_SNR40(t) = ((conc_DPFadd5_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
            perc_error_BB_DPFadd15_Eadd15_SNR40(t) = ((conc_DPFadd15_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFadd10_Eadd10_SNR40(t) = ((conc_DPFadd10_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFadd5_Eadd5_SNR40(t) = ((conc_DPFadd5_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFtake15_Etake15_SNR40(t) = ((conc_DPFtake15_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFtake10_Etake10_SNR40(t) = ((conc_DPFtake10_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFtake5_Etake5_SNR40(t) = ((conc_DPFtake5_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFtake15_Eadd10_SNR40(t) = ((conc_DPFtake15_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFtake15_Eadd5_SNR40(t) = ((conc_DPFtake15_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFtake10_Eadd15_SNR40(t) = ((conc_DPFtake10_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFtake10_Eadd5_SNR40(t) = ((conc_DPFtake10_Eadd5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFadd5_Etake15_SNR40(t) = ((conc_DPFadd5_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFadd5_Etake10_SNR40(t) = ((conc_DPFadd5_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFtake5_Eadd15_SNR40(t) = ((conc_DPFtake5_Eadd15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFtake5_Eadd10_SNR40(t) = ((conc_DPFtake5_Eadd10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFadd10_Etake15_SNR40(t) = ((conc_DPFadd10_Etake15_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFadd10_Etake5_SNR40(t) = ((conc_DPFadd10_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_DPFadd15_Etake10_SNR40(t) = ((conc_DPFadd15_Etake10_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            perc_error_BB_DPFadd15_Etake5_SNR40(t) = ((conc_DPFadd15_Etake5_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            
            perc_error_BB_5s_SNR40(t) = max(abs([perc_error_BB_DPFtake5_Etake5_SNR40(t), perc_error_BB_DPFadd5_Eadd5_SNR40(t), perc_error_BB_DPFadd5_Etake5_SNR40(t), perc_error_BB_DPFtake5_Eadd5_SNR40(t)]));
            perc_error_BB_10s_SNR40(t) = max(abs([perc_error_BB_DPFtake10_Etake10_SNR40(t), perc_error_BB_DPFadd10_Eadd10_SNR40(t), perc_error_BB_DPFadd10_Etake10_SNR40(t), perc_error_BB_DPFtake10_Eadd10_SNR40(t)]));
            perc_error_BB_15s_SNR40(t) = max(abs([perc_error_BB_DPFtake15_Etake15_SNR40(t), perc_error_BB_DPFadd15_Eadd15_SNR40(t), perc_error_BB_DPFadd15_Etake15_SNR40(t), perc_error_BB_DPFtake15_Eadd15_SNR40(t)]));
        end
        %%
        perc_accum_errors_DPFE_BB = [perc_diff_BB; perc_error_BB_5s_SNR40; perc_error_BB_10s_SNR40; perc_error_BB_15s_SNR40];
        
        % FigH = figure('Position', get(0, 'Screensize'));
        % F = getframe(FigH);
        % imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient and DPF percentage changes for broadband system 1.png', 'png')
        % plot(percent_vector,perc_accum_errors_DPFE_BB(:,1),'r*-','MarkerSize',20)
        % hold on
        % plot(percent_vector,perc_accum_errors_DPFE_BB(:,2),'bs-','MarkerSize',10)
        % plot(percent_vector,perc_accum_errors_DPFE_BB(:,5),'go-','MarkerSize',10)
        % ax = gca;
        % ax.FontSize = 20;
        % %ylim([0 35]);
        % legend('HbO','HbR','oxCCO','Location','Best');
        % xlabel('Percentage change in both extinction coefficient and DPF values');
        % ylabel('Percentage error in chromophore concentration values');
        % title({'Comparison of percentage errors for increasing extinction coefficient', 'and DPF percentage changes for broadband system 1 (680-921nm)'});
        % saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1.png')
        
        
        
        %
        %
        %
        %---------------------------------------------------------------------------------
        %% Code for LED
        % x contains the selected wavelengths
        % x = [874,889,892,895,898];
        x_setup = linspace(700,900,nwavs);
        x = ceil(x_setup);
        
        % Bandwidth
        w_LED = zeros(length(x),1) + 0.1;    
        % Broadband spectrum
        wl = 680:921; %nm
        % wl = 771:906; %nm
        
        % Load chromophore data
        % Generate vectors containing extinction coefficient and mua data
        % Using Journal of Biomedical Optics 15 5, 056002 September/October 2010
        % GetExtinctions(lambda)
        % Returns the specific absorption coefficients for [HbO Hb H2O lipid aa3]
        % for the specified wavelengths. Note that the specific
        % absorption coefficient (defined base e) is equal to the 
        % specific extinction coefficient (defined base 10) times 2.303.
        E = GetExtinctions(x); % (5 wavs)
        E_new = GetExtinctions(wl); % (broadband)
        
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
        
        E_rewrite = [E(:,1) E(:,2)]; % ADD
        E = E_rewrite;
        
        E_new_rewrite = [E_new(:,1) E_new(:,2)]; % ADD
        E_new = E_new_rewrite;
        
        %% Lorentzian
        
        % Approximate LED shape as a Lorentzian function using different bandwidths
        for i = 1:length(wl)
            for j = 1:length(x)
                LED_spec(i,j) = (1/pi)*((w_LED(j)/2)/(((wl(i)-x(j))^2)+(w_LED(j)/2)^2));
            end
        end
        
        % %% Chromophore spectra
        % %Plot constituent spectra
        % figure;
        % plot(x,E(:,1),'ro-');hold on;
        % plot(x,E(:,2),'bo-');
        % plot(x,(E(:,3)/10),'go-'); %Divide by 10 to improve display
        % plot(x,(E(:,4)/10),'ko-'); %Divide by 10 to improve display
        % plot(x,E(:,5),'mo-');
        % xlim([min(x) max(x)])
        % legend('HbO','HbR','H2O','Lipid','aa3');
        % xlabel('Wavelength (nm)');
        % ylabel('Specific abs. coeff or Abs. coeff/10 (mm^{-1}M^{-1} or mm^{-1})');
        % title('Constituent spectra (and polynomial fits)');
        
        % Find absorption and intensity for broadband case
        % muabrain_broad_1 = (E_new(:,1)*C_HbO(1)) + (E_new(:,2)*C_HbR(1)) + (E_new(:,3)*W(1)) + (E_new(:,4)*L1(1)) + (E_new(:,5)*C_aa3(1));
        muabrain_broad_1 = (E_new(:,1)*C_HbO(1)) + (E_new(:,2)*C_HbR(1));
        
        I_broadband_1 = exp(-muabrain_broad_1);
        
        % Find area under each LED peak
        LED_area = trapz(LED_spec);
        
        % For LED
        for l = 1:length(x)
            for k = 1:length(wl)
                I_effective_1_LED(k,l) = change_I_BB(k)*LED_spec(k,l);
                %I_effective_1_LED(k,l) = I_broadband_1(k)*LED_spec(k,l);
            end
            I_effective_noise_LED(:,l) = awgn(I_effective_1_LED(:,l),60); %???
        end
        
        for g = 1:length(x)
            I_prime_1_LED(g) = sum(I_effective_1_LED(:,g))/LED_area(g);
            %I_prime_1_LED(g) = sum(I_effective_noise_LED(:,g))/LED_area(g);
        end
        
        %% Effective absorption
        
        % Find the extinction coefficients adjusted by convolution according to LED
        % bandwidth
        % For LED
        for l = 1:length(ac_conc_change) %chromophores 
            for j = 1:length(x) % wavelengths
                for k = 1:length(wl) % broadband
                    e_effective(k,l) = E_new(k,l)*LED_spec(k,j); % column = chromophore, row = wavelength
                end
                E_prime_LED(j,l) = sum(e_effective(:,l))/LED_area(l);
            end
        end
        
        E_prime_invert_LED = pinv(E_prime_LED);
        
        % % Again for extinction vals from book
        % for l = 1:length(x) %chromophores 
        %     for j = 1:length(x) % wavelengths
        %         for k = 1:length(wl) % broadband
        %             e_effective_book(k,l) = ext_BioBook_new(k,l)*LED_spec(k,j); % column = chromophore, row = wavelength
        %         end
        %         E_prime_LED_book(j,l) = sum(e_effective_book(:,l))/LED_area(l);
        %     end
        % end
        % 
        % E_prime_invert_LED_book = pinv(E_prime_LED_book);
        
        % Find effective wavelengths used when LED shifted by Lorentzian
        % (wavelength value associated with new extinction coefficient)
        tolerance = 10; % +-5nm
        
        % For LED
        
        for u = 1:length(concs1) 
            for r = 1:length(x)
                [val, idx] = min(abs(E_new((x(r)-wl(1))-tolerance:(x(r)-wl(1))+tolerance,u) - E_prime_LED(r,u)));
                new_ext_wavs_init(r,u) = idx;
                new_ext_wavs_LED(r,u) = new_ext_wavs_init(r,u) + (x(r)-tolerance);
            end
        end
        
        for r = 1:length(x)
            ave_ext_wav_LED(:,r) = sum(new_ext_wavs_LED(r,:))/length(concs1);
        end
        
        ave_ext_wav_LED = round(ave_ext_wav_LED);
        
        for g = 1:length(ave_ext_wav_LED)
            new_ave_DPF_LED(g) = ave_DPF_BB(ave_ext_wav_LED(g)-wl(1));
        end
        
        % % Exts from book
        % 
        % for u = 1:length(concs1) 
        %     for r = 1:length(x)
        %         [val, idx] = min(abs(E_new((x(r)-wl(1))-tolerance:(x(r)-wl(1))+tolerance,u) - E_prime_invert_LED_book(r,u)));
        %         new_ext_wavs_init(r,u) = idx;
        %         new_ext_wavs_book_LED(r,u) = new_ext_wavs_init(r,u) + (x(r)-tolerance);
        %     end
        % end
        % 
        % for r = 1:length(x)
        %     ave_ext_wav_book_LED(:,r) = sum(new_ext_wavs_book_LED(r,:))/length(x);
        % end
        % 
        % ave_ext_wav_book_LED = round(ave_ext_wav_book_LED);
        % 
        % for g = 1:length(ave_ext_wav_LED)
        %     new_ave_DPF_book_LED(g) = ave_DPF_BB(ave_ext_wav_book_LED(g)-wl(1));
        % end
        
        % Get input mua for diffusion code according to LED spectra
        % mua1_LED = (E_prime_LED(:,1)*C_HbO(1)) + (E_prime_LED(:,2)*C_HbR(1)) + (E_prime_LED(:,3)*W(1)) + (E_prime_LED(:,4)*L1(1)) + (E_prime_LED(:,5)*C_aa3(1));
        mua1_LED = (E_prime_LED(:,1)*C_HbO(1)) + (E_prime_LED(:,2)*C_HbR(1));
        
        % This is second input to diffusion equations
        % mua2_LED = (E_prime_LED(:,1)*C_HbO(2)) + (E_prime_LED(:,2)*C_HbR(2)) + (E_prime_LED(:,3)*W(2)) + (E_prime_LED(:,4)*L1(2)) + (E_prime_LED(:,5)*C_aa3(2));
        mua2_LED = (E_prime_LED(:,1)*C_HbO(2)) + (E_prime_LED(:,2)*C_HbR(2));
        
            %% Slab model from Contini et al. Appl.Opt. 36, 4587 (1997). - 1
        %     
            power = -1.2;
            musp = 1.0 * (x/800).^power;
        %     figure
        %     plot(x,musp,'ko-');
        %     xlim([min(x) max(x)])
        %     xlabel('Wavelength (nm)');
        %     ylabel('Transport scatter coefficient (mm^{-1})');
        %     title('Transport scatter coefficient LED');
            
            mut = mua1_LED + musp';
            D = 1 ./ (3 * mut');
            mueff = sqrt(3 .* mua1_LED' .* mut');
            z0 = 1 ./ musp;
        
            % Mismatched boundary parameters
            nrel = 1.4; % Refractive index of medium
            Apar = 504.332889-(2641.00214*nrel)+(5923.699064*nrel^2)-(7376.335814*nrel^3)+(5507.53041*nrel^4)-(2463.357945*nrel^5)+(610.956547*nrel^6)-(64.8047*nrel^7);
            ze = (2 * Apar) .* D;
        
            % Sum terms to generate reflectance: Equation (45)
            mlimit = 4;
            R1_LED = zeros([1 length(x)]); 
            %R(l) = num2cell(R);
            for m =(-mlimit):1:(mlimit)
                z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
                z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
                g = (SD * SD) + (z3m .* z3m);
                h = (SD * SD) + (z3m .* z3m);
                rmg = sqrt(mua1_LED' .* g ./ D);
                rmh = sqrt(mua1_LED' .* h ./ D);
                ermg = exp(-rmg);
                ermh = exp(-rmh);
                RT11 = z3m .* (g.^-1.5) .* (rmg + 1) .* ermg;
                RT21 = z4m .* (h.^-1.5) .* (rmh + 1) .* ermh;
                R1_LED = R1_LED + ((RT11 - RT21) ./ (-4 * pi));
        %         R{l} = R;
            end
                %R{l} = num2cell(R);
            %R1(1,:) = R1;
            % Sum terms to generate DPF: Equation (47)
            mlimit = 4;
            L1_LED = zeros([1 length(x)]); 
            for m =(-mlimit):1:(mlimit)
                z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
                z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
                g = (SD * SD) + (z3m .* z3m);
                h = (SD * SD) + (z3m .* z3m);
                rmg = sqrt(mua1_LED' .* g ./ D);
                rmh = sqrt(mua1_LED' .* h ./ D);
                ermg = exp(-rmg);
                ermh = exp(-rmh);
                LT0 = (-8 * pi) .* D .* R1_LED;
                LT11 = z3m .* (g.^-0.5) .* ermg;
                LT21 = z4m .* (h.^-0.5) .* ermh;
                L1_LED = L1_LED + ((LT11 - LT21) ./ LT0); 
            end
            %L1_LED = awgn(L1_LED,5,'measured');
          
        
        %     %% Calculate and plot the DPF
            DPF1_LED = L1_LED ./ SD; % Wavelength dependent
        %     figure
        %     plot(x,DPF1_LED,'ko-');
        %     xlim([min(x) max(x)])
        %     xlabel('Wavelength (nm)');
        %     ylabel('DPF');
        %     title('DPF 1 LED');
        % 
        %     %% Plot attenuation (-log diffuse reflectance)
            A1 = -log(R1_LED);
        %     figure
        %     subplot(1,2,1);
        %     plot(x,R1_LED,'bo-');
        %     xlim([min(x) max(x)]);
        %     xlabel('Wavelength (nm)');
        %     ylabel('Diffuse Reflectance');
        %     title('Diffuse Reflectance 1 LED');
        %     subplot(1,2,2);
        %     plot(x,A1,'ko-');
        %     xlim([min(x) max(x)]);
        %     xlabel('Wavelength (nm)');
        %     ylabel('-log_{e}Reflectance');
        %     title('Attenuation = -ln(Reflectance) 1 LED');
        
            %% Slab model from Contini et al. Appl.Opt. 36, 4587 (1997). - 2
            
            mut = mua2_LED + musp';
            D = 1 ./ (3 * mut');
            mueff = sqrt(3 .* mua2_LED' .* mut');
            z0 = 1 ./ musp;
        
            % Mismatched boundary parameters
            nrel = 1.4; % Refractive index of medium
            Apar = 504.332889-(2641.00214*nrel)+(5923.699064*nrel^2)-(7376.335814*nrel^3)+(5507.53041*nrel^4)-(2463.357945*nrel^5)+(610.956547*nrel^6)-(64.8047*nrel^7);
            ze = (2 * Apar) .* D;
        
            % Sum terms to generate reflectance: Equation (45)
            mlimit = 4;
            R2_LED = zeros([1 length(x)]); 
            %R(l) = num2cell(R);
            for m =(-mlimit):1:(mlimit)
                z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
                z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
                g = (SD * SD) + (z3m .* z3m);
                h = (SD * SD) + (z3m .* z3m);
                rmg = sqrt(mua2_LED' .* g ./ D);
                rmh = sqrt(mua2_LED' .* h ./ D);
                ermg = exp(-rmg);
                ermh = exp(-rmh);
                RT12 = z3m .* (g.^-1.5) .* (rmg + 1) .* ermg;
                RT22 = z4m .* (h.^-1.5) .* (rmh + 1) .* ermh;
                R2_LED = R2_LED + ((RT12 - RT22) ./ (-4 * pi));
        %         R{l} = R;
            end
                %R{l} = num2cell(R);
            % Sum terms to generate DPF: Equation (47)
            mlimit = 4;
            L2_LED = zeros([1 length(x)]); 
            for m =(-mlimit):1:(mlimit)
                z3m = - z0 - (4 .* m .* ze) - (2 .* m .* slab);
                z4m = z0 - (((4 .* m) - 2) .* ze) - (2 .* m .* slab);  
                g = (SD * SD) + (z3m .* z3m);
                h = (SD * SD) + (z3m .* z3m);
                rmg = sqrt(mua2_LED' .* g ./ D);
                rmh = sqrt(mua2_LED' .* h ./ D);
                ermg = exp(-rmg);
                ermh = exp(-rmh);
                LT0 = (-8 * pi) .* D .* R2_LED;
                LT12 = z3m .* (g.^-0.5) .* ermg;
                LT22 = z4m .* (h.^-0.5) .* ermh;
                L2_LED = L2_LED + ((LT12 - LT22) ./ LT0); 
            end
            
        
        %     %% Calculate and plot the DPF
            DPF2_LED = L2_LED ./ SD; % Wavelength dependent
        %     figure
        %     plot(x,DPF2_LED,'ko-');
        %     xlim([min(x) max(x)])
        %     xlabel('Wavelength (nm)');
        %     ylabel('DPF');
        %     title('DPF 2 LED');
        % 
        %     %% Plot attenuation (-log diffuse reflectance)
            A2 = -log(R2_LED);
        %     figure
        %     subplot(1,2,1);
        %     plot(x,R2_LED,'bo-');
        %     xlim([min(x) max(x)]);
        %     xlabel('Wavelength (nm)');
        %     ylabel('Diffuse Reflectance');
        %     title('Diffuse Reflectance 2 LED');
        %     subplot(1,2,2);
        %     plot(x,A2,'ko-');
        %     xlim([min(x) max(x)]);
        %     xlabel('Wavelength (nm)');
        %     ylabel('-log_{e}Reflectance');
        %     title('Attenuation = -ln(Reflectance) 2 LED');
        
        %% Task 2
        
        change_I_LED = log(R2_LED) - log(R1_LED);
        
        SNR_vector = [20; 30; 40; 50; 60];
        
        for l = 1:20
        
            R1_LED_SNR20 = awgn(R1_LED,20,'measured');
            R2_LED_SNR20 = awgn(R2_LED,20,'measured');
            change_I_LED_SNR20 = log(R2_LED_SNR20) - log(R1_LED_SNR20);
            % change_I_LED_SNR20 = log(R2_LED_SNR20) - log(R1_LED);
        
            R1_LED_SNR30 = awgn(R1_LED,30,'measured');
            R2_LED_SNR30 = awgn(R2_LED,30,'measured');
            change_I_LED_SNR30 = log(R2_LED_SNR30) - log(R1_LED_SNR30);
            % change_I_LED_SNR30 = log(R2_LED_SNR30) - log(R1_LED);
        
            R1_LED_SNR40 = awgn(R1_LED,40,'measured');
            R2_LED_SNR40 = awgn(R2_LED,40,'measured');
            change_I_LED_SNR40 = log(R2_LED_SNR40) - log(R1_LED_SNR40);
            % change_I_LED_SNR40 = log(R2_LED_SNR40) - log(R1_LED);
        
            R1_LED_SNR50 = awgn(R1_LED,50,'measured');
            R2_LED_SNR50 = awgn(R2_LED,50,'measured');
            change_I_LED_SNR50 = log(R2_LED_SNR50) - log(R1_LED_SNR50);
            % change_I_LED_SNR50 = log(R2_LED_SNR50) - log(R1_LED);
        
            R1_LED_SNR60 = awgn(R1_LED,60,'measured');
            R2_LED_SNR60 = awgn(R2_LED,60,'measured');
            change_I_LED_SNR60 = log(R2_LED_SNR60) - log(R1_LED_SNR60);
            % change_I_LED_SNR60 = log(R2_LED_SNR60) - log(R1_LED);
        
            ave_DPF_LED = (DPF1_LED + DPF2_LED)/2;
        
            for t = 1:length(x)
                calc_I_LED(t) = ave_DPF_LED(t)*SD*(C_HbO(2)-C_HbO(1))*E_prime_LED(t,1);
            end
        
            for t = 1:length(x)
                perc_diff_I_LED(t) = 100 - ((calc_I_LED(t)/change_I_LED(t))*100);
            end
        
            % This is the value to pay attention to!! From diffusion
            conc_LED =  E_prime_invert_LED * (-change_I_LED ./ (SD*ave_DPF_LED))'; % From diffusion
        
            conc_LED_SNR20 =  E_prime_invert_LED * (-change_I_LED_SNR20 ./ (SD*ave_DPF_LED))';
            conc_LED_SNR30 =  E_prime_invert_LED * (-change_I_LED_SNR30 ./ (SD*ave_DPF_LED))';
            conc_LED_SNR40 =  E_prime_invert_LED * (-change_I_LED_SNR40 ./ (SD*ave_DPF_LED))';
            conc_LED_SNR50 =  E_prime_invert_LED * (-change_I_LED_SNR50 ./ (SD*ave_DPF_LED))';
            conc_LED_SNR60 =  E_prime_invert_LED * (-change_I_LED_SNR60 ./ (SD*ave_DPF_LED))';
        
            for t = 1:length(ac_conc_change)
                perc_diff_LED(t) = ((conc_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                perc_diff_LED_SNR20(t) = ((conc_LED_SNR20(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                perc_diff_LED_SNR30(t) = ((conc_LED_SNR30(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                perc_diff_LED_SNR40(t) = ((conc_LED_SNR40(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                perc_diff_LED_SNR50(t) = ((conc_LED_SNR50(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                perc_diff_LED_SNR60(t) = ((conc_LED_SNR60(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
            end
        
            perc_vector_SNR = [20;30;40;50;60;100];
            perc_diff_SNR_vector = [perc_diff_LED_SNR20;perc_diff_LED_SNR30;perc_diff_LED_SNR40;perc_diff_LED_SNR50;perc_diff_LED_SNR60;perc_diff_LED];
            perc_diff_SNR_vector = abs(perc_diff_SNR_vector);
            
            if l == 1 && pert_number == 1
        
%  
                %% Introduce percentage errors in LED measurements
        
                % Introduce a 5% error in DPF
                percent_vector = [0.01; 5; 10; 15];
        
                for r = 1:length(ave_DPF_LED)
                    five_perc_LED_DPF(r) = (ave_DPF_LED(r)/100) *5;
                    ten_perc_LED_DPF(r) = (ave_DPF_LED(r)/100) *10;
                    fifteen_perc_LED_DPF(r) = (ave_DPF_LED(r)/100) *15;
        
                    LED_DPF_add5(r) = ave_DPF_LED(r) + five_perc_LED_DPF(r);
                    LED_DPF_add10(r) = ave_DPF_LED(r) + ten_perc_LED_DPF(r);
                    LED_DPF_add15(r) = ave_DPF_LED(r) + fifteen_perc_LED_DPF(r);
        
                    LED_DPF_take5(r) = ave_DPF_LED(r) - five_perc_LED_DPF(r);
                    LED_DPF_take10(r) = ave_DPF_LED(r) - ten_perc_LED_DPF(r);
                    LED_DPF_take15(r) = ave_DPF_LED(r) - fifteen_perc_LED_DPF(r);
                end
        
                conc_LED_add5percDPF =  E_prime_invert_LED * (-change_I_LED ./ (SD*LED_DPF_add5))';
                conc_LED_add10percDPF =  E_prime_invert_LED * (-change_I_LED ./ (SD*LED_DPF_add10))';
                conc_LED_add15percDPF =  E_prime_invert_LED * (-change_I_LED ./ (SD*LED_DPF_add15))';
        
                conc_LED_take5percDPF =  E_prime_invert_LED * (-change_I_LED ./ (SD*LED_DPF_take5))';
                conc_LED_take10percDPF =  E_prime_invert_LED * (-change_I_LED ./ (SD*LED_DPF_take10))';
                conc_LED_take15percDPF =  E_prime_invert_LED * (-change_I_LED ./ (SD*LED_DPF_take15))';
        
                for t = 1:length(ac_conc_change)
                    perc_error_LED_add5DPF(t) = ((conc_LED_add5percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_add10DPF(t) = ((conc_LED_add10percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_add15DPF(t) = ((conc_LED_add15percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_take5DPF(t) = ((conc_LED_take5percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_take10DPF(t) = ((conc_LED_take10percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_take15DPF(t) = ((conc_LED_take15percDPF(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    % take largest value
                    perc_error_LED_DPF5(t) = max(perc_error_LED_add5DPF(t),perc_error_LED_take5DPF(t));
                    perc_error_LED_DPF10(t) = max(perc_error_LED_add10DPF(t),perc_error_LED_take10DPF(t));
                    perc_error_LED_DPF15(t) = max(perc_error_LED_add15DPF(t),perc_error_LED_take15DPF(t));
                end
        
                perc_error_DPF_vector_LED = [perc_diff_LED; perc_error_LED_DPF5; perc_error_LED_DPF10; perc_error_LED_DPF15];
                perc_error_DPF_vector_LED = abs(perc_error_DPF_vector_LED);
        
    %             FigH = figure('Position', get(0, 'Screensize'));
    %             F = getframe(FigH);
    %             imwrite(F.cdata, 'Comparison of percentage errors for increasing DPF percentage changes for 5-wavelength LED system.png', 'png')
    %             plot(percent_vector,perc_error_DPF_vector_LED(:,1),'r*-','MarkerSize',20)
    %             hold on
    %             plot(percent_vector,perc_error_DPF_vector_LED(:,2),'bs-','MarkerSize',10)
    %             % plot(percent_vector,perc_error_DPF_vector_LED(:,3),'go-','MarkerSize',10)
    %             ax = gca;
    %             ax.FontSize = 20;
    %             ylim([0 35]);
    %             legend('HbO','HbR','Location','Best');
    %             xlabel('Percentage change in DPF value');
    %             ylabel('Percentage error in chromophore concentration values');
    %             title({'Comparison of percentage errors for increasing', 'DPF percentage changes for 5-wavelength LED system'});
    %             saveas(gcf,'Comparison of percentage errors for increasing DPF percentage changes for 5-wavelength LED system.png')
        
                for r = 1:length(E_prime_LED)
                    for q = 1:length(ac_conc_change)
                        five_perc_LED_E(r,q) = (E_prime_LED(r,q)/100) *5;
                        ten_perc_LED_E(r,q) = (E_prime_LED(r,q)/100) *10;
                        fifteen_perc_LED_E(r,q) = (E_prime_LED(r,q)/100) *15;
        
                        LED_E_add5(r,q) = E_prime_LED(r,q) + five_perc_LED_E(r,q);
                        LED_E_add10(r,q) = E_prime_LED(r,q) + ten_perc_LED_E(r,q);
                        LED_E_add15(r,q) = E_prime_LED(r,q) + fifteen_perc_LED_E(r,q);
        
                        LED_E_take5(r,q) = E_prime_LED(r,q) - five_perc_LED_E(r,q);
                        LED_E_take10(r,q) = E_prime_LED(r,q) - ten_perc_LED_E(r,q);
                        LED_E_take15(r,q) = E_prime_LED(r,q) - fifteen_perc_LED_E(r,q);
                    end
                end
        
                conc_LED_add5percE =  pinv(LED_E_add5) * (-change_I_LED ./ (SD*ave_DPF_LED))';
                conc_LED_add10percE =  pinv(LED_E_add10) * (-change_I_LED ./ (SD*ave_DPF_LED))';
                conc_LED_add15percE =  pinv(LED_E_add15) * (-change_I_LED ./ (SD*ave_DPF_LED))';
        
                conc_LED_take5percE =  pinv(LED_E_take5) * (-change_I_LED ./ (SD*ave_DPF_LED))';
                conc_LED_take10percE =  pinv(LED_E_take10) * (-change_I_LED ./ (SD*ave_DPF_LED))';
                conc_LED_take15percE =  pinv(LED_E_take15) * (-change_I_LED ./ (SD*ave_DPF_LED))';
        
                for t = 1:length(ac_conc_change)
                    perc_error_LED_add5E(t) = ((conc_LED_add5percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_add10E(t) = ((conc_LED_add10percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_add15E(t) = ((conc_LED_add15percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_take5E(t) = ((conc_LED_take5percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_take10E(t) = ((conc_LED_take10percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_take15E(t) = ((conc_LED_take15percE(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_E5(t) = max(perc_error_LED_add5E(t),perc_error_LED_take5E(t));
                    perc_error_LED_E10(t) = max(perc_error_LED_add10E(t),perc_error_LED_take10E(t));
                    perc_error_LED_E15(t) = max(perc_error_LED_add15E(t),perc_error_LED_take15E(t));
                end
        
                perc_error_E_vector_LED = [perc_diff_LED; perc_error_LED_E5; perc_error_LED_E10; perc_error_LED_E15];
                perc_error_E_vector_LED = abs(perc_error_E_vector_LED);
        
    %             FigH = figure('Position', get(0, 'Screensize'));
    %             F = getframe(FigH);
    %             imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient percentage changes for 5-wavelength LED system.png', 'png')
    %             plot(percent_vector,perc_error_E_vector_LED(:,1),'r*-','MarkerSize',20)
    %             hold on
    %             plot(percent_vector,perc_error_E_vector_LED(:,2),'bs-','MarkerSize',10)
    %             % plot(percent_vector,perc_error_E_vector_LED(:,3),'go-','MarkerSize',10)
    %             ax = gca;
    %             ax.FontSize = 20;
    %             ylim([0 35]);
    %             legend('HbO','HbR','Location','Best');
    %             xlabel('Percentage change in extinction coefficient value');
    %             ylabel('Percentage error in chromophore concentration values');
    %             title({'Comparison of percentage errors for increasing extinction coefficient', 'percentage changes for 5-wavelength LED system'});
    %             saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for 5-wavelength LED system.png')
            
            end
        
                %% Again for SNR 50
        
                %
                conc_DPFadd5_Eadd5_SNR50_LED = pinv(LED_E_add5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add5))';
        
                conc_DPFadd10_Eadd5_SNR50_LED = pinv(LED_E_add5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add10))';
        
                conc_DPFadd15_Eadd5_SNR50_LED = pinv(LED_E_add5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add15))';
        
                %
                conc_DPFadd10_Eadd10_SNR50_LED = pinv(LED_E_add10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add10))';
        
                conc_DPFadd15_Eadd15_SNR50_LED = pinv(LED_E_add15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add15))';
        
                %
                %------------------------------------------------------------
                %
                conc_DPFtake5_Etake5_SNR50_LED = pinv(LED_E_take5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take5))';
        
                conc_DPFtake10_Etake5_SNR50_LED = pinv(LED_E_take5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take10))';
        
                conc_DPFtake15_Etake5_SNR50_LED = pinv(LED_E_take5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take15))';
        
                %
                conc_DPFtake10_Etake10_SNR50_LED = pinv(LED_E_take10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take10))';
        
                conc_DPFtake15_Etake15_SNR50_LED = pinv(LED_E_take15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take15))';
        
                %
                %--------------------------------------------------------------------
                %
                conc_DPFadd5_Etake5_SNR50_LED = pinv(LED_E_take5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add5))';
        
                conc_DPFadd10_Etake5_SNR50_LED = pinv(LED_E_take5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add10))';
        
                conc_DPFadd15_Etake5_SNR50_LED = pinv(LED_E_take5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add15))';
        
                %
                conc_DPFadd10_Etake10_SNR50_LED = pinv(LED_E_take10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add10))';
        
                conc_DPFadd15_Etake15_SNR50_LED = pinv(LED_E_take15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add15))';
        
                %-------------------------------------------------------------
                %
                conc_DPFtake5_Eadd5_SNR50_LED = pinv(LED_E_add5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take5))';
        
                conc_DPFtake10_Eadd5_SNR50_LED = pinv(LED_E_add5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take10))';
        
                conc_DPFtake15_Eadd5_SNR50_LED = pinv(LED_E_add5) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take15))';
        
                %
                conc_DPFtake10_Eadd10_SNR50_LED = pinv(LED_E_add10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take10))';
        
                conc_DPFtake15_Eadd15_SNR50_LED = pinv(LED_E_add15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take15))';
        
                conc_DPFtake15_Eadd10_SNR50_LED = pinv(LED_E_add10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take15))';
        
                conc_DPFtake10_Eadd15_SNR50_LED = pinv(LED_E_add15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take10))';
        
                conc_DPFadd5_Etake15_SNR50_LED = pinv(LED_E_take15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add5))';
        
                conc_DPFadd5_Etake10_SNR50_LED = pinv(LED_E_take10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add5))';
        
                conc_DPFtake5_Eadd15_SNR50_LED = pinv(LED_E_add15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take5))';
        
                conc_DPFtake5_Eadd10_SNR50_LED = pinv(LED_E_add10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_take5))';
        
                conc_DPFadd10_Etake15_SNR50_LED = pinv(LED_E_take15) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add10))';
        
                conc_DPFadd15_Etake10_SNR50_LED = pinv(LED_E_take10) * (-change_I_LED_SNR50 ./ (SD*LED_DPF_add15))';
        
                for t = 1:length(ac_conc_change)
                    perc_error_LED_DPFtake15_Eadd15_SNR50(t) = ((conc_DPFtake15_Eadd15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFtake10_Eadd10_SNR50(t) = ((conc_DPFtake10_Eadd10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFtake5_Eadd5_SNR50(t) = ((conc_DPFtake5_Eadd5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFadd15_Etake15_SNR50(t) = ((conc_DPFadd15_Etake15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFadd10_Etake10_SNR50(t) = ((conc_DPFadd10_Etake10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFadd5_Etake5_SNR50(t) = ((conc_DPFadd5_Etake5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFadd15_Eadd15_SNR50(t) = ((conc_DPFadd15_Eadd15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFadd10_Eadd10_SNR50(t) = ((conc_DPFadd10_Eadd10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFadd5_Eadd5_SNR50(t) = ((conc_DPFadd5_Eadd5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFtake15_Etake15_SNR50(t) = ((conc_DPFtake15_Etake15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFtake10_Etake10_SNR50(t) = ((conc_DPFtake10_Etake10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFtake5_Etake5_SNR50(t) = ((conc_DPFtake5_Etake5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFtake15_Eadd10_SNR50(t) = ((conc_DPFtake15_Eadd10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFtake15_Eadd5_SNR50(t) = ((conc_DPFtake15_Eadd5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFtake10_Eadd15_SNR50(t) = ((conc_DPFtake10_Eadd15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFtake10_Eadd5_SNR50(t) = ((conc_DPFtake10_Eadd5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFadd5_Etake15_SNR50(t) = ((conc_DPFadd5_Etake15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFadd5_Etake10_SNR50(t) = ((conc_DPFadd5_Etake10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFtake5_Eadd15_SNR50(t) = ((conc_DPFtake5_Eadd15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFtake5_Eadd10_SNR50(t) = ((conc_DPFtake5_Eadd10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFadd10_Etake15_SNR50(t) = ((conc_DPFadd10_Etake15_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFadd10_Etake5_SNR50(t) = ((conc_DPFadd10_Etake5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_DPFadd15_Etake10_SNR50(t) = ((conc_DPFadd15_Etake10_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
                    perc_error_LED_DPFadd15_Etake5_SNR50(t) = ((conc_DPFadd15_Etake5_SNR50_LED(t) - ac_conc_change(t))/ac_conc_change(t)) * 100;
        
                    perc_error_LED_5s_SNR50(t) = max(abs([perc_error_LED_DPFtake5_Etake5_SNR50(t), perc_error_LED_DPFadd5_Eadd5_SNR50(t), perc_error_LED_DPFadd5_Etake5_SNR50(t), perc_error_LED_DPFtake5_Eadd5_SNR50(t)]));
                    perc_error_LED_10s_SNR50(t) = max(abs([perc_error_LED_DPFtake10_Etake10_SNR50(t), perc_error_LED_DPFadd10_Eadd10_SNR50(t), perc_error_LED_DPFadd10_Etake10_SNR50(t), perc_error_LED_DPFtake10_Eadd10_SNR50(t)]));
                    perc_error_LED_15s_SNR50(t) = max(abs([perc_error_LED_DPFtake15_Etake15_SNR50(t), perc_error_LED_DPFadd15_Eadd15_SNR50(t), perc_error_LED_DPFadd15_Etake15_SNR50(t), perc_error_LED_DPFtake15_Eadd15_SNR50(t)]));
                end
                %%
                perc_accum_errors_DPFE_LEDSNR50 = [perc_diff_LED; perc_error_LED_5s_SNR50; perc_error_LED_10s_SNR50; perc_error_LED_15s_SNR50];
                perc_accum_errors_DPFE_LEDSNR50 = abs(perc_accum_errors_DPFE_LEDSNR50);
    
                perc_accum_errors_DPFE_LEDSNR50_saved{l} = perc_accum_errors_DPFE_LEDSNR50(:,:); %Save each iteration to cell array
        end
                % FigH = figure('Position', get(0, 'Screensize'));
                % F = getframe(FigH);
                % imwrite(F.cdata, 'Comparison of percentage errors for increasing extinction coefficient and DPF percentage changes for broadband system 1 SNR 50.png', 'png')
                % plot(percent_vector,perc_accum_errors_DPFE_LEDSNR50(:,1),'r*-','MarkerSize',20)
                % hold on
                % plot(percent_vector,perc_accum_errors_DPFE_LEDSNR50(:,2),'bs-','MarkerSize',10)
                % plot(percent_vector,perc_accum_errors_DPFE_LEDSNR50(:,5),'go-','MarkerSize',10)
                % ax = gca;
                % ax.FontSize = 20;
                % %ylim([0 35]);
                % legend('HbO','HbR','oxCCO','Location','Best');
                % xlabel('Percentage change in both extinction coefficient and DPF values');
                % ylabel('Percentage error in chromophore concentration values');
                % title({'Comparison of percentage errors for increasing extinction coefficient', 'and DPF percentage changes for broadband system 1 (680-921nm), SNR = 50'});
                % saveas(gcf,'Comparison of percentage errors for increasing extinction coefficient percentage changes for broadband system 1 SNR 50.png')
            for u = 1:4
                for j = 1:2
                    for k = 1:10
                        vector_error_LEDSNR50(k,j) = perc_accum_errors_DPFE_LEDSNR50_saved{1,k}(u,j);
                    end
                    ave_vector_error_LEDSNR50(u,j) = mean(vector_error_LEDSNR50(:,j));
                    high_bar_vector_error_LEDSNR50(u,j) = max(vector_error_LEDSNR50(:,j));
                    low_bar_vector_error_LEDSNR50(u,j) = min(vector_error_LEDSNR50(:,j));
                    stand_dev_LEDSNR50(u,j) = std(vector_error_LEDSNR50(:,j));
                end
            end
    
            pert_ave_error_LED{pert_number} = ave_vector_error_LEDSNR50(:,:);
            pert_high_bar_error_LED{pert_number} = high_bar_vector_error_LEDSNR50;
            pert_low_bar_error_LED{pert_number} = low_bar_vector_error_LEDSNR50;
            pert_stand_dev_LED{pert_number} = stand_dev_LEDSNR50;
    
    end
        
        %%
        
        %%
        % Average LED
        for u = 1:4
            for j = 1:2
                for k = 1:pert_number
                    total_pert_error_LEDSNR50(k,j) = pert_ave_error_LED{1,k}(u,j);
                    high_bar_vector_error_LEDSNR50(k,j) = pert_high_bar_error_LED{1,k}(u,j);
                    low_bar_vector_error_LEDSNR50(k,j) = pert_low_bar_error_LED{1,k}(u,j);
                    standard_dev_LEDSNR50(k,j) = pert_stand_dev_LED{1,k}(u,j);
                end
                val_to_ave_error = total_pert_error_LEDSNR50(:,j);
                high_to_ave_error = high_bar_vector_error_LEDSNR50(:,j);
                low_to_ave_error = low_bar_vector_error_LEDSNR50(:,j);
                stdtoave = standard_dev_LEDSNR50(:,j);
                ave_pert_error_LEDSNR50(u,j) = mean(val_to_ave_error(val_to_ave_error<Inf));
                high_bar_pert_error_LEDSNR50(u,j) = mean(high_to_ave_error(high_to_ave_error<Inf));
                low_bar_pert_error_LEDSNR50(u,j) = mean(low_to_ave_error(low_to_ave_error<Inf));
                final_std_val_LEDSNR50(u,j) = mean(stdtoave(stdtoave<Inf));
            end
        end
    
        num_wavs = length(x);
        save(['ave_error_',num2str(num_wavs),'_wavelengths_LUMO'],'ave_pert_error_LEDSNR50','final_std_val_LEDSNR50')
        
        clear all;
    
end
    %%
% % Graph for LUMO
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Average percentage errors for increasing E and DPF changes for LUMO SNR 50 random conc.png', 'png')
% errorbar(percent_vector,ave_pert_error_LEDSNR50(:,1),low_bar_pert_error_LEDSNR50(:,1),high_bar_pert_error_LEDSNR50(:,1),'r*-','MarkerSize',20,'LineWidth',1.8)
% hold on
% errorbar(percent_vector,ave_pert_error_LEDSNR50(:,2),low_bar_pert_error_LEDSNR50(:,2),high_bar_pert_error_LEDSNR50(:,2),'bs-','MarkerSize',20,'LineWidth',1.5)
% % errorbar(percent_vector,ave_pert_error_LEDSNR50(:,3),low_bar_pert_error_LEDSNR50(:,3),high_bar_pert_error_LEDSNR50(:,3),'go-','MarkerSize',20,'LineWidth',1)
% grid on
% ax = gca;
% ax.FontSize = 20;
% xlim([-5 20]);
% legend('HbO','HbR','Location','Best');
% xlabel('Percentage change in both extinction coefficient and DPF values');
% ylabel('Average percentage error in chromophore concentration values');
% title({'Average percentage errors for increasing extinction coefficient', 'and DPF percentage changes for LUMO with random conc changes, SNR = 50'});
% saveas(gcf,'Average percentage errors for increasing E and DPF changes for LUMO SNR 50 random conc.png')

