
% Enter simulation data
SD = 30; % Source-detector separation on same surface of slab (mm)
slab = 50; % Slab thickness (mm)

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

E_rewrite = [E(:,1) E(:,2) E(:,5)]; % ADD
E = E_rewrite;

E_new_rewrite = [E_new(:,1) E_new(:,2) E_new(:,5)]; % ADD
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


% Plot
% FigH = figure('Position', get(0, 'Screensize'));
% F = getframe(FigH);
% imwrite(F.cdata, 'Constituent spectra (and polynomial fits) from Homer2 for broadband system 1.png', 'png')
plot(x,E(:,1),'ro-');hold on;
plot(x,E(:,2),'bo-');
% plot(x,(E(:,3)/10),'go-'); %Divide by 10 to improve display
% plot(x,(E(:,4)/10),'ko-'); %Divide by 10 to improve display
plot(x,E(:,3),'mo-');
xlim([min(x) max(x)])
legend('HbO','HbR','aa3');
xlabel('Wavelength (nm)');
ylabel('Specific abs. coeff or Abs. coeff/10 (mm^{-1}M^{-1} or mm^{-1})');
title('Constituent spectra from Homer2 for broadband system 1 (680-921nm)');
% saveas(gcf,'Constituent spectra (and polynomial fits) from Homer2 for broadband system 1.png')


