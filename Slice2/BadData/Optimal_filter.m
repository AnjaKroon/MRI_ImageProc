%% Creates K-Space data after WH Optimal Filtering %%

clc; clear all; close all;

%%%%%% LOAD BAD DATA %%%%%%
load('slice2_channel1_bad.mat');
load('slice2_channel2_bad.mat');
load('slice2_channel3_bad.mat');

%%%%%% LOAD GOOD DATA %%%%%%
load('/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/GoodData/slice2_channel1_good.mat');
load('/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/GoodData/slice2_channel2_good.mat');
load('/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/GoodData/slice2_channel3_good.mat');
whos

%%%%%% RUN OPT FILT %%%%%%
slice2_channel1_optFilt = Op_filter(slice2_channel1_goodData,slice2_channel1_badData,10,64);
slice2_channel2_optFilt = Op_filter(slice2_channel2_goodData,slice2_channel2_badData,10,64);
slice2_channel3_optFilt = Op_filter(slice2_channel3_goodData,slice2_channel3_badData,10,64);

%%%%%% CHECK OBTAINED RESULTS %%%%%%
imagesc(100*log(abs(slice2_channel1_optFilt)));
imagesc(100*log(abs(slice2_channel2_optFilt)));
imagesc(100*log(abs(slice2_channel3_optFilt)));

%%%%%% SAVE OBTAINED RESULTS %%%%%%
save('slice2_channel1_optFilt.mat', 'slice2_channel1_optFilt');
save('slice2_channel2_optFilt.mat', 'slice2_channel2_optFilt');
save('slice2_channel3_optFilt.mat', 'slice2_channel3_optFilt');

%%%%%% OPTIMUM FILTERING CODE %%%%%%
function X = Op_filter(G,B,O,p)               % enter good K-space image, bad image ,filter order and k-space patch size
    dimx= size(G,1);
    dimy= size(G,2);
    Rg = Autocorrelation(G,dimx,dimy,p);      % calculate the autocorrelation values ot the K-subspaces
    Rb = Autocorrelation(B,dimx,dimy,p);

    for i = 1: size(Rg,1)
        for j =1: size(Rg,2)
            g = cell2mat(Rg(i,j));
            b = cell2mat(Rb(i,j));
            g = g(1:O);      
            b = b(1:O);
            R_x = toeplitz(b');
            W{i,j} = inv(R_x)*g';             % calculate the optimal filter coefficients for each K-subspace
        end
    end
    x= zeros(dimx/p,1)+p;                     % define the x dimensions of the k-subspaces
    y= zeros(dimy/p,1)+p;                     % define the y dimensions of the k-subspaces
    B= mat2cell(B,x,y);

    for k = 1:size(W,1)
        for l= 1:size(W,2)
            for m = 1:p
                w= cell2mat(W(k,l));
                b= cell2mat(B(k,l));
                f=conv(b(m,:),w);                  % filter each row of the K-subspace of the bad image with optimal coefficients
                F(m,:)= f(1:p);
                Fc{k,l} = F;
           end
        end
    end
    X=cell2mat(Fc);

end

















