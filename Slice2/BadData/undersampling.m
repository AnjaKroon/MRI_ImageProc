clc; clear all; close all;

load('slice2_channel1_bad.mat');
load('slice2_channel2_bad.mat');
load('slice2_channel3_bad.mat');
whos

%%%%%%%%%%%%%%% BERNOULLI DENSITY %%%%%%%%%%%%%%%

% how do you know what p value to set?

p = 0.5; 
berMask = rand(size(slice2_channel1_badData)) < p;

slice2_channel1_berKS = slice2_channel1_badData .* berMask;
slice2_channel2_berKS = slice2_channel2_badData .* berMask;
slice2_channel3_berKS = slice2_channel3_badData .* berMask;

%%%%%%%%%%%%%%% GAUSSIAN DENSITY %%%%%%%%%%%%%%%
% variable density sampling
% more concentrated towards the middle

% how do you know what mean and var to set

sigma = 1; 
mean = 0; 

numCols = size(slice2_channel1_badData, 2);
numRows = size(slice2_channel1_badData, 1);

stepsizeX = 2 / (numCols - 1);
stepsizeY = 2 / (numRows - 1);

xvalues = -1:stepsizeX:1;
yvalues = -1:stepsizeY:1;

[x, y] = meshgrid(xvalues, yvalues);

gauMask = exp(-(x-mean).^2 / (2*sigma^2)) .* exp(-(y-mean).^2 / (2*sigma^2));

gauMask = gauMask / max(gauMask(:)); % normalize

slice2_channel1_gauKS = slice2_channel1_badData .* gauMask;
slice2_channel2_gauKS = slice2_channel2_badData .* gauMask;
slice2_channel3_gauKS = slice2_channel3_badData .* gauMask;



%%%%%%%%%%%%%%% EXPONENTIAL DENSITY %%%%%%%%%%%%%%%
% how do you know what alpha value to set
alpha = 0.1; 

numCols = size(slice2_channel1_badData, 2);
numRows = size(slice2_channel1_badData, 1);

stepsizeX = 2 / (numCols - 1);
stepsizeY = 2 / (numRows - 1);

xvalues = -1:stepsizeX:1;
yvalues = -1:stepsizeY:1;

[x, y] = meshgrid(xvalues, yvalues);

expMask = exp(-alpha * sqrt(x.^2 + y.^2));

slice2_channel1_expKS = slice2_channel1_badData .* expMask;
slice2_channel2_expKS = slice2_channel2_badData .* expMask;
slice2_channel3_expKS = slice2_channel3_badData .* expMask;

%%%%%%%%%%%%%%% CARTESIAN GRID %%%%%%%%%%%%%%%
lines = 2;

numCols = size(slice2_channel1_badData, 2);
numRows = size(slice2_channel1_badData, 1);

carMask = zeros(numRows, numCols);
inx = 1:lines:numRows;
carMask(inx, :) = 1;

slice2_channel1_carKS = slice2_channel1_badData .* carMask;
slice2_channel2_carKS = slice2_channel2_badData .* carMask;
slice2_channel3_carKS = slice2_channel3_badData .* carMask;


berMask = berMask(1:128, 1:128);
carMask = carMask(1:128, 1:128);

%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%
figure;

subplot(1, 2, 1);
imshow(berMask, []);
title('Ber Mask');
colormap(gray);
colorbar;

subplot(1, 2, 2);
imshow(carMask, []);
title('Car Mask');
colormap(gray);
colorbar;

axis off;

%%%%%%%%%%%%%%% SAVING KSPACE DATA %%%%%%%%%%%%%%%


filts = {'ber', 'gau', 'exp', 'car'};
disp(size(filts))
for i=1:length(filts)
    for j = 1:3
        matfileName = strcat('slice2_channel',string(j),"_", string(filts{i}), 'KS.mat');
        varName = strcat('slice2_channel',string(j),"_", string(filts{i}), 'KS');
        %disp(matfileName)
        %disp(varName)
        save(matfileName, varName);
    end
end

