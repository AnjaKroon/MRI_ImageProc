%% Comparing wRSS and RSS approaches using evaluation metrics %%

clear all; close all;

%%%%%%% PULL IN GOOD DATA %%%%%%%
load('slice2_image_goodData.mat');
whos

%%%%%%% PULL IN MODIFIED DATA %%%%%%%
badPath = '/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/BadData';
badContents = dir(fullfile(badPath, 'slice2_spatial*'));
spatialKData = cell(1, min(3, length(badContents)));

%%%%%%% LOAD DATA %%%%%%%
for i = 1:length(spatialKData)
    filePath = fullfile(badPath, badContents(i).name);
    loadedData = load(filePath);  
    spatialKData{i} = loadedData
end

%spatialKData{1} % regular
%spatialKData{2} % weighted

spatial_reg = cell2mat(struct2cell(spatialKData{1}))
spatial_w = cell2mat(struct2cell(spatialKData{2}))

load('slice2_image_goodData.mat');

% perform SNR comparison
psnr_reg = psnr(outputImage, spatial_reg)
psnr_w = psnr(outputImage, spatial_w)

% perform SSIM
ssim_reg = ssim(outputImage, spatial_reg)
ssim_w = ssim(outputImage, spatial_w)

% perform MSE
nmse_reg = nmse(outputImage, spatial_reg)
nmse_w = nmse(outputImage, spatial_w)

function nmseValue = nmse(orig, comp)
    mse = sum((orig - comp).^2) / length(orig);
    origEnergy = sum(orig.^2);
    nmseValue = mse / origEnergy;
end