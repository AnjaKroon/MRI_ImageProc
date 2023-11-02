%% Compares restore/fuse method vs. fuse/restore method %%

clc; clear all; close all;

%%%%%% LOAD BAD DATA %%%%%%
load('slice2_channel1_bad.mat');
load('slice2_channel2_bad.mat');
load('slice2_channel3_bad.mat');

%% Weighted RSS %%

%%%%%%%%% ESTIMATING NOISE FROM EACH CHANNEL %%%%%%%%%
eye_ch_1(:,:) = (ifftshift(ifft2(slice2_channel1_badData),1));
eye_ch_2(:,:) = (ifftshift(ifft2(slice2_channel2_badData),1));
eye_ch_3(:,:) = (ifftshift(ifft2(slice2_channel3_badData),1));

not_cropped_ch1 = eye_ch_1;
not_cropped_ch2 = eye_ch_2;
not_cropped_ch3 = eye_ch_3;

crop_x = [128 + 60 : 348 - 33]; 
eye_ch_1 = eye_ch_1(crop_x, :);
eye_ch_2 = eye_ch_2(crop_x, :);
eye_ch_3 = eye_ch_3(crop_x, :);

eye_ch_1 = reshape(squeeze(eye_ch_1(:,:)),[128 128]); 
eye_ch_2 = reshape(squeeze(eye_ch_2(:,:)),[128 128]); 
eye_ch_3 = reshape(squeeze(eye_ch_3(:,:)),[128 128]); 

black_ch_1 = eye_ch_1;
black_ch_2 = eye_ch_2;
black_ch_3 = eye_ch_3;


%%%%%%% CH1: APPLY INTENSITY ADJUSTMENTS %%%%%%%
std_within = 0.995; 
% set maximum intensity to contain 99.5 % of intensity values per image
[aa, val] = hist(eye_ch_1(:),linspace(0,max(...
                                    eye_ch_1(:)),1000));
    thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
% set threshold value to 65536
eye_ch_1 = uint16(eye_ch_1 * 65536 / thresh);

black_ch_1 = black_ch_1([1:10], [119:128]);
black_ch_1 = double(black_ch_1 * 65536 / thresh);


%%%%%%% CH2: APPLY INTENSITY ADJUSTMENTS %%%%%%%
[aa, val] = hist(eye_ch_2(:),linspace(0,max(...
                                    eye_ch_2(:)),1000));
    thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
% set threshold value to 65536

eye_ch_2 = uint16(eye_ch_2 * 65536 / thresh);

black_ch_2 = black_ch_2([1:10], [119:128]);
black_ch_2 = double(black_ch_2 * 65536 / thresh);


%%%%%%% CH3: APPLY INTENSITY ADJUSTMENTS %%%%%%%
[aa, val] = hist(eye_ch_3(:),linspace(0,max(...
                                    eye_ch_3(:)),1000));
    thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
% set threshold value to 65536
eye_ch_3 = uint16(eye_ch_3 * 65536 / thresh);

black_ch_3 = black_ch_3([1:10], [119:128]);
black_ch_3 = double(black_ch_3 * 65536 / thresh);


%%%%%%% CREATING GROUND TRUTH FOR BLACK SQUARE %%%%%%%
ref = ones(10);
ref = double(ref);

err_ch1 = immse(black_ch_1, ref);
err_ch2 = immse(black_ch_2, ref);
err_ch3 = immse(black_ch_3, ref);

var_ch1 = var(black_ch_1)
var_ch2 = var(black_ch_2)
var_ch3 = var(black_ch_3)

noise_ch1 = black_ch_1 - ref; % noise from ch 1
noise_ch2 = black_ch_2 - ref;
noise_ch3 = black_ch_3 - ref;

psnr_ch1 = psnr(black_ch_1, ref)
psnr_ch2 = psnr(black_ch_2, ref)
psnr_ch3 = psnr(black_ch_3, ref)

norm_psnr_factor = min([(psnr_ch1), (psnr_ch2), (psnr_ch3)]);

weights = [psnr_ch1/norm_psnr_factor, psnr_ch2/norm_psnr_factor, psnr_ch3/norm_psnr_factor]


% clear compensation, preparation, based on fourier transformed blinked 
% k-space data (Data_raw)
clear_comp = linspace(10,0.1,size(eye_ch_1,2)).^2; 
clear_matrix = repmat(clear_comp,[size(eye_ch_1,1) 1]);

% combine 3 channels sum of squares and add clear compensation
weighted_sos  = sqrt( abs(squeeze(weights(1) .* not_cropped_ch1(:,:,1))).^2 + ...
            abs(weights(2) .* squeeze(not_cropped_ch2(:,:,1))).^2 + ...
           abs(squeeze(weights(3) .* not_cropped_ch3(:,:,1))).^2).* clear_matrix; 

% crop images because we are only interested in eye. Make it square 
% 128 x 128
crop_x = [128 + 60 : 348 - 33]; % crop coordinates
weighted_sos = weighted_sos(crop_x, :);

% Visualize the images. 
%image
eye_weighted = reshape(squeeze(weighted_sos(:,:)),[128 128]); 
% For better visualization and contrast of the eye images, histogram based
% compensation will be done 
std_within = 0.995; 
% set maximum intensity to contain 99.5 % of intensity values per image
[aa, val] = hist(eye_weighted(:),linspace(0,max(...
                                    eye_weighted(:)),1000));
    thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
% set threshold value to 65536
eye_weighted = uint16(eye_weighted * 65536 / thresh);

slice2_spatial_SoSw = eye_weighted;
save("slice2_spatial_SoSw.mat", "slice2_spatial_SoSw")


%% Regular RSS %%

% IFFT of k-space data
%channel 1 (replace "slice1_channel1_goodData" with
%slice1_channel1_badData) for bad images
Data_img(:,:,1) = ifftshift(ifft2(slice2_channel1_badData),1);
%channel 2
Data_img(:,:,2) = ifftshift(ifft2(slice2_channel2_badData),1);
%channel 3
Data_img(:,:,3) = ifftshift(ifft2(slice2_channel2_badData),1);

% clear compensation, preparation, based on fourier transformed blinked 
% k-space data (Data_raw)
clear_comp = linspace(10,0.1,size(Data_img,2)).^2; 
clear_matrix = repmat(clear_comp,[size(Data_img,1) 1]);

% combine 3 channels sum of squares and add clear compensation
eye_raw  = sqrt( abs(squeeze(Data_img(:,:,1))).^2 + ...
           abs(squeeze(Data_img(:,:,2))).^2 + ...
           abs(squeeze(Data_img(:,:,3))).^2).* clear_matrix; 

% crop images because we are only interested in eye. Make it square 
% 128 x 128
crop_x = [128 + 60 : 348 - 33]; % crop coordinates
eye_raw = eye_raw(crop_x, :);

% Visualize the images. 
%image
eye_visualize_orig = reshape(squeeze(eye_raw(:,:)),[128 128]); 
% For better visualization and contrast of the eye images, histogram based
% compensation will be done 
std_within = 0.995; 
% set maximum intensity to contain 99.5 % of intensity values per image
[aa, val] = hist(eye_visualize_orig(:),linspace(0,max(...
                                    eye_visualize_orig(:)),1000));
    thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
% set threshold value to 65536
eye_visualize_orig = uint16(eye_visualize_orig * 65536 / thresh);

slice2_spatial_SoS = eye_visualize_orig;
save("slice2_spatial_SoS.mat", "slice2_spatial_SoS")



%% plotting scripts
%close all
%figure(1); 
%imagesc(eye_ch_1(:,:));
%axis image, 
%colormap gray;
%axis off

%figure(2); 
%imagesc(eye_ch_2(:,:));
%axis image, 
%colormap gray;
%axis off

%figure(3); 
%imagesc(eye_ch_3(:,:));
%axis image, 
%colormap gray;
%axis off

figure(4); 
imagesc(eye_weighted(:,:));
title("Weighted SoS")
axis image, 
colormap gray;
axis off

figure(5); 
imagesc(eye_visualize_orig(:,:));
title("Regular SoS")
axis image, 
colormap gray;
axis off

