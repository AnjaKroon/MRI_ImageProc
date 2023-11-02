clc; clear all; close all;

%%%%%%% LOADING DATA %%%%%%%
load('slice2_channel1_bad.mat');
load('slice2_channel2_bad.mat');
load('slice2_channel3_bad.mat');
whos

% README: Replace with path of local files
load('/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/GoodData/slice2_channel1_good.mat');
load('/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/GoodData/slice2_channel2_good.mat');
load('/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/GoodData/slice2_channel3_good.mat');
whos

%%%%%%% CALLING CUBIC INTERPOLATION %%%%%%%
slice2_channel1_cub = Q_int(slice2_channel1_badData)
slice2_channel2_cub = Q_int(slice2_channel2_badData)
slice2_channel3_cub = Q_int(slice2_channel3_badData)

% imagesc(100*log(abs(slice2_channel1_cub)));


%%%%%%% SAVING LOCALLY %%%%%%%
save('slice2_channel1_cub.mat', 'slice2_channel1_cub');
save('slice2_channel2_cub.mat', 'slice2_channel2_cub');
save('slice2_channel3_cub.mat', 'slice2_channel3_cub');

function I = Q_int(G)  
    M = size(G,2);              % input K-space image
    N = size(G,1);
    [T,x,y] = SL_detection(G,M); % detect the corrupted spectral lines
    y(y== 65) = [];
    disp(x)
    G(:, y) = 0;                 % setting the corrupted columns to zero 
   
    for i = 1:length(y)                                   
        for j = 1:N
           B = [(y(i)-2)^3 (y(i)-2)^2 y(i)-2 1; (y(i)-1)^3 (y(i)-1)^2 y(i)-1 1; (y(i)+1)^3 (y(i)+1)^2 y(i)+1 1; (y(i)+2)^3 (y(i)+2)^2 y(i)+2 1 ];
           F = [G(j,y(i)-2); G(j,y(i)-1); G(j,y(i)+1); G(j,y(i)+2)];
           inv(B);
           A = inv(B)*F;
           X = [y(i)^3 y(i)^2 y(i) 1];
           G(j,y(i))= A.'*X.';                          % replacing the corrupted columns by interpolated versions
        end
    end
    I=G;

    
end

