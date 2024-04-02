%% Setup initial values

theFiles = dir('/Users/georginaleadley/Documents/GitHub/Diffusion-Simulations/Code replica/Official_50000_combs/*.mat');

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
    data_to_analyse = newdata.master_list_sorted;
    
    best_vals{k} = data_to_analyse(1,:);
    middle_val = length(data_to_analyse)/2;
    middle_int = ceil(middle_val);
    mid_vals{k} = data_to_analyse(middle_int,:);
    worst_val = length(data_to_analyse)*0.95;
    worst_int = ceil(worst_val);
    worst_vals{k} = data_to_analyse(worst_int,:);
end

save(['cond_array_3-50wav'],'best_vals','mid_vals','worst_vals')

