clc; clear all; close all;

load('slice2_channel1.mat');
load('slice2_channel2.mat');
load('slice2_channel3.mat');

% for temporary testing purposes
load('slice2_channel1_gauKS.mat')
load('slice2_channel1_gauKS.mat')
load('slice2_channel1_gauKS.mat')
whos

% B - k-space bad
% G - k-space good
% M - 512
    % N - 128
    % O - order
    % p - 32

o = 10;
p = 32;

slice2_channel1_pred = Pf(slice2_channel1_badData, slice2_channel1_gauKS, 512, 128, o, p)
%slice2_channel2_pred = Pf(slice2_channel2_badData, slice2_channel1_gauKS, 512, 128, o, p)
%slice2_channel3_pred = Pf(slice2_channel3_badData, slice2_channel1_gauKS, 512, 128, o, p)

imagesc(100*log(abs(slice2_channel1_pred)));

function Prediction_filter= Pf(B,G,M,N,O,p)
    % B - k-space bad
    % G - k-space good
    % M - 512
    % N - 128
    % O - order
    % p - 32

    [T,x,y] = SL_detection(B,N);    % obtain the corrupted spectral lines 
    disp(y)
    y(y == 65) = [];                                                                    
    R_b = Autocorrelation(B,M,N,p);  % compute the autocorrelation values of the K-subspaces
    R_g = Autocorrelation(G,M,N,p)   

    for i = 1:length(y)      % identify to which subspace each corruptes value corresponds
        for j = 1:M
            if(floor(M/p)~=M/p)
            cell_row = floor(M/p)+1;
            else 
            cell_row = floor(M/p);
            end
            if(floor(y(i)/p)~=y(i)/p)
            cell_column = floor(y(i)/p)+1;
            else 
            cell_column = floor(y(i)/p);
            end
            b=cell2mat(R_b(cell_row,cell_column)); % obtain the required autocorrelation values
            g=cell2mat(R_g(cell_row,cell_column));

            y1= b(1:O);
            R_x = toeplitz(y1');
            R_a= g(2:O+1);  
            W = inv(R_x)*R_a';  % compute the filter coefficients
            for k= 1:O
            A(k) = W(k)*B(j,y(i)+k);                                        
            end
            B(j,y(i))= sum(A); % replace the corrupted values by the predicted ones
        end
    end
    Prediction_filter = B;

    function [T,x,y] = SL_detection(K, M) % input corrupted K-space image
        % M= size(K,2);
        K= sqrt(abs(K));   
                                                                
        for i = 1:M-2  % compute the sum of the difference between adjacent columns 
           T(i)= sum(abs(K(:,i+1)-K(:,i) + K(:,i+1)-K(:,i+2))); 
        end
    
        T = [0, T];
        [x,y] = findpeaks(T,'MinPeakProminence',75); % detect large differences between adjacent columns 
    end
end






