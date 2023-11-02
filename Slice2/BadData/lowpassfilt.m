%% Creates the filtered k space data for all LPFs %%

% filter implementations
    
    clear all;
    clc; clear all; close all;

    %%%%%% LOAD BAD DATA %%%%%%
    load('slice2_channel1_bad.mat');
    load('slice2_channel2_bad.mat');
    load('slice2_channel3_bad.mat');
    whos
   
    %%%%%% CREATE FILTERS %%%%%%
    gauLPF = gaussian_lpf(120,128,512);

    buterLPF = butterworth_lpf(4,120,128,512);

    idealLPF = ideal_lpf(120, 128, 512);


    %%%%%% CONFIRM FILTER VISUALLY %%%%%%
    figure;
    imshow(buterLPF, []);
    title('Butterworth LPF (Order 4, Sigma 120)');
    colormap(gray);
    colorbar;

    %%%%%% APPLY FILTER %%%%%%

    slice2_channel1_gauLPF = gauLPF .* slice2_channel1_badData;
    slice2_channel2_gauLPF = gauLPF .* slice2_channel2_badData;
    slice2_channel3_gauLPF = gauLPF .* slice2_channel3_badData;

    slice2_channel1_butterLPF = buterLPF .* slice2_channel1_badData;
    slice2_channel2_butterLPF = buterLPF .* slice2_channel2_badData;
    slice2_channel3_butterLPF = buterLPF .* slice2_channel3_badData;

    slice2_channel1_idealLPF = idealLPF .* slice2_channel1_badData;
    slice2_channel2_idealLPF = idealLPF .* slice2_channel2_badData;
    slice2_channel3_idealLPF = idealLPF .* slice2_channel3_badData;

    %%%%%% SAVE FILTERED DATA %%%%%%
    save('slice2_channel1_gauLPF.mat', 'slice2_channel1_gauLPF');
    save('slice2_channel2_gauLPF.mat', 'slice2_channel2_gauLPF');
    save('slice2_channel3_gauLPF.mat', 'slice2_channel3_gauLPF');

    save('slice2_channel1_butterLPF.mat', 'slice2_channel1_butterLPF');
    save('slice2_channel2_butterLPF.mat', 'slice2_channel2_butterLPF');
    save('slice2_channel3_butterLPF.mat', 'slice2_channel3_butterLPF');

    save('slice2_channel1_idealLPF.mat', 'slice2_channel1_idealLPF');
    save('slice2_channel2_idealLPF.mat', 'slice2_channel2_idealLPF');
    save('slice2_channel3_idealLPF.mat', 'slice2_channel3_idealLPF');
    

    %%%%%% FUNCTIONS %%%%%%
    function H = butterworth_lpf(n,F,M,N) % to approach the     .
        u=0:(N-1);
        v=0:(M-1);
        idx=find(u>N/2);
        u(idx)=u(idx)-N;
        idy=find(v>M/2);
        v(idy)=v(idy)-M;
        [V,U]=meshgrid(v,u);
        D=sqrt(U.^2+V.^2);
        H = fftshift(1./(1 + (D./F).^(2*n)));
        %imshow((H),[])
    end
    function H = gaussian_lpf(F,M,N) % no ringing due to smoothness IDFT is also gaussian
        u=0:(N-1);
        v=0:(M-1);
        idx=find(u>N/2);
        u(idx)=u(idx)-N;
        idy=find(v>M/2);
        v(idy)=v(idy)-M;
        [V,U]=meshgrid(v,u);
        D=sqrt(U.^2+V.^2);
        H = fftshift(exp(-(D.^2)./(2*(F^2))))
        %imshow((H),[])
    end      
    function H = ideal_lpf(F,M,N) % unfortunally ringing
        u=0:(N-1);
        v=0:(M-1);
        idx=find(u>N/2);
        u(idx)=u(idx)-N;
        idy=find(v>M/2);
        v(idy)=v(idy)-M;
        [V,U]=meshgrid(v,u);
        D=sqrt(U.^2+V.^2);
        H = fftshift(double(D<=F));
        %imshow((H),[])
    end