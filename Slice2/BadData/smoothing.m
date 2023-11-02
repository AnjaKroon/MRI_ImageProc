clc; clear all; close all;

%%%%%% LOADING BAD DATA %%%%%%
load('slice2_channel1_bad.mat');
load('slice2_channel2_bad.mat');
load('slice2_channel3_bad.mat');

%%%%%%%%%%%%%%%%%%%%% PERFORM SMOOTHING %%%%%%%%%%%%%%%%%%%%%

%%%%%%%% Linear Method 1: Gaussian filter %%%%%%%%
sigma = 1;   % larger sigma, more smoothing
% badData = {slice2_channel1_badData, slice2_channel2_badData, slice2_channel3_badData};

slice2_channel1_gauss = [];
slice2_channel2_gauss = [];
slice2_channel3_gauss = [];

% Perform gaussian filtering by separating the real and im. components
for i = 1:3
    if i == 1
        realPart = real(slice2_channel1_badData);
        imaginaryPart = imag(slice2_channel1_badData);
        smoothedReal = imgaussfilt(realPart, sigma);
        smoothedImaginary = imgaussfilt(imaginaryPart, sigma);
        slice2_channel1_gauss = complex(smoothedReal, smoothedImaginary);
        %save('slice2_allch_gauss.mat', 'slice2_channel1_gauss');
    elseif i == 2
        realPart = real(slice2_channel2_badData);
        imaginaryPart = imag(slice2_channel2_badData);
        smoothedReal = imgaussfilt(realPart, sigma);
        smoothedImaginary = imgaussfilt(imaginaryPart, sigma);
        slice2_channel2_gauss = complex(smoothedReal, smoothedImaginary);
        %save('slice2_allch_gauss.mat', 'slice2_channel2_gauss');
    else
        realPart = real(slice2_channel3_badData);
        imaginaryPart = imag(slice2_channel3_badData);
        smoothedReal = imgaussfilt(realPart, sigma);
        smoothedImaginary = imgaussfilt(imaginaryPart, sigma);
        slice2_channel3_gauss = complex(smoothedReal, smoothedImaginary);
        %save('slice2_allch_gauss.mat', 'slice2_channel3_gauss');
    end
end


save('slice2_channel1_gauss.mat', 'slice2_channel1_gauss');
save('slice2_channel2_gauss.mat', "slice2_channel2_gauss");
save('slice2_channel3_gauss.mat', "slice2_channel3_gauss");

% to check if saving works properly
% load('slice2_allch_gauss.mat')
% whos

% to check if smoothing is doing something
% diff = slice2_channel1_modData - slice2_channel1_badData;
% disp(diff)

% to check if filtering works
%close all
%figure(1); 
%imagesc(100*log(abs(slice2_channel1_badData)));
%title('Origional Image')
%figure(2); 
%imagesc(100*log(abs(slice2_channel1_modData_gauss)));
%title('Gaussian Filter Smoothed Image');


%%%%%%%%% Linear Method 2: Mean filter %%%%%%%%

% Mean filter size (3x3, 5x5, etc.)
filterSize = 3;

% Perform mean filter smoothing, see function at end of file
slice2_channel1_mean = meanFilter(slice2_channel1_badData, ...
    filterSize);
slice2_channel2_mean = meanFilter(slice2_channel2_badData, ...
    filterSize);
slice2_channel3_mean = meanFilter(slice2_channel3_badData, ...
    filterSize);

% to check: display one channel with filter
% figure(3); 
% imagesc(100*log(abs(slice2_channel1_badData)));
% title("Original Image")
% figure(4); 
% imagesc(100*log(abs(slice2_channel1_mean)));
% title('Mean Filter Smoothed Image');

% save three variables in slice2_allch_mean.mat
save('slice2_channel1_mean.mat', 'slice2_channel1_mean');
save('slice2_channel2_mean.mat', "slice2_channel2_mean");
save('slice2_channel3_mean.mat', "slice2_channel3_mean");



%%%%%%%%% Non-Linear Method 1: Median filter %%%%%%%%
% Mean filter size (3x3, 5x5, etc.)
filterSize = 3;

% Perform mean filter smoothing, see function at end of file
slice2_channel1_med = medianFilter(slice2_channel1_badData, ...
    filterSize);
slice2_channel2_med = medianFilter(slice2_channel2_badData, ...
    filterSize);
slice2_channel3_med = medianFilter(slice2_channel3_badData, ...
    filterSize);

% to check: display one channel with filter
% to check: display one channel with filter
figure(5); 
imagesc(100*log(abs(slice2_channel1_badData)));
title("Original Image")
figure(6); 
imagesc(100*log(abs(slice2_channel1_med)));
title('Median Filter Smoothed Image');

% to check: display one channel with filter
figure(5); 
imagesc(100*log(abs(slice2_channel1_badData)));
title("Original Image")
figure(6); 
imagesc(100*log(abs(slice2_channel1_med)));
title('Median Filter Smoothed Image');

figure(5); 
imagesc(100*log(abs(slice2_channel1_badData)));
title("Original Image")
figure(6); 
imagesc(100*log(abs(slice2_channel1_med)));
title('Median Filter Smoothed Image');

% save three variables in slice2_allch_mean.mat
save('slice2_channel1_med.mat', 'slice2_channel1_med');
save('slice2_channel2_med.mat', "slice2_channel2_med");
save('slice2_channel3_med.mat', "slice2_channel3_med");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean filter function
function output = meanFilter(input, filterSize)
    paddedInput = padarray(input, [filterSize, filterSize], 0, 'both');
    [rows, cols] = size(paddedInput);
    halfFilterSize = floor(filterSize / 2);
    output = zeros(rows, cols);
    for i = 1 + halfFilterSize : rows - halfFilterSize
        for j = 1 + halfFilterSize : cols - halfFilterSize
            % Extract local neighborhood
            neighborhood = paddedInput(i - halfFilterSize : i + ...
                halfFilterSize, j - halfFilterSize : j + halfFilterSize);
            % Apply mean filter
            output(i, j) = mean(neighborhood(:));
        end
    end
    output = output(filterSize:end-1-filterSize, filterSize:end-1-filterSize);
end


% Median filter function
function output = medianFilter(input, filterSize)
    paddedInput = padarray(input, [filterSize, filterSize], 0, 'both');
    [rows, cols] = size(paddedInput);
    halfFilterSize = floor(filterSize / 2);
    output = zeros(rows, cols);
    for i = 1 + halfFilterSize : rows - halfFilterSize
        for j = 1 + halfFilterSize : cols - halfFilterSize
            % Extract local neighborhood
            neighborhood = paddedInput(i - halfFilterSize : i + ...
                halfFilterSize, j - halfFilterSize : j + halfFilterSize);
            % Apply mean filter
            output(i, j) = median(neighborhood(:));
        end
    end
    output = output(filterSize:end-1-filterSize, filterSize:end-1-filterSize);
end



