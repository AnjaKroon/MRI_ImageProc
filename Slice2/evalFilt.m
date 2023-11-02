function evaluateFiltering = evalFilt(good, filt)
    
    %%%%%%%%%%%%%%%% IMPORT GOOD K DATA %%%%%%%%%%%%%%%%
    % Pull in three good k space data files 
    goodPath = '/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/GoodData';
    goodContents = dir(fullfile(goodPath, '*.mat'))
    goodKData = cell(1, min(3, length(goodContents)));  

    for i = 1:length(goodKData)
        filePath = fullfile(goodPath, goodContents(i).name);
        loadedData = load(filePath);  
        goodKData{i} = loadedData;
    end
    
    %%%%%%%%%%%%%%%% IMPORT ALL FILTERED K DATA %%%%%%%%%%%%%%%%
    % Pull in filtered k space data files
    badPath = '/Users/anja/Documents/Delft/Q1/SDSPM/project/MRI_datasets/Slice2/BadData';
    badContents = dir(fullfile(badPath, 'slice2_channel*'))
    filtKData = cell(1, min(46, length(badContents))); 

    % List of methods used. Important for their order
    methods = {'bad', 'berKS', 'butterLPF', 'carKS', 'cub', 'expKS', 'gauKS', 'gauLPF', 'gauss', 'idealLPF', 'mean', 'med', 'optFilt', 'pred'}

    for i = 1:length(filtKData)
        % disp(badContents(i).name);
        filePath = fullfile(badPath, badContents(i).name);
        loadedData = load(filePath);
        filtKData{i} = loadedData;
    end
    %disp(size(filtKData))
   

    %%%%%%%%%%%%%%%% CREATE ALL EYE IMAGES %%%%%%%%%%%%%%%%
    % generate image from good k space data
    % and for all filters
        % generate image from the output of filtered k space
    
    % turn structs into matricies
    goodMat = {[], [], []};
    for i = 1:length(goodKData)
        cellArray = struct2cell(goodKData{i});
        goodMat{i} = cell2mat(cellArray); % index 1 is ch 1 etc. 
    end
    
    % turn good k space data into an image
    good_eye_test = getImage(goodMat{1}, goodMat{2}, goodMat{3});
    outputImage = good_eye_test(:,:,1);
    save("slice2_image_goodData.mat", 'outputImage');
    
    % prepare an empty array to be populated with the filtered k space data
    rows = length(methods);
    cols = 3;
    filtKSMat = cell(rows, cols);
    for i = 1:rows
        for j = 1:cols
            filtKSMat{i, j} = [];
        end
    end
    
    n = rows % number of rows in filtMat
    for i = 1:length(filtKData) %42
        % bad1, bad2, bad3 etc.
        %col = mod(i,4)
        %if col == 0
            %floor(i/n)
            %col = n;
        %end
        %row_ = floor((i-0.1)/n)+1
        cellArray = struct2cell(filtKData{i});
        filtKData{i} = cell2mat(cellArray);  
    end
    size(filtKData)
    filtKSMat = reshape(filtKData, n+1, 3);
    % disp(size(filtKSMat))

    % disp(filtKSMat)

    % for each row in filtMat, generate the corresponding image, and save
    for i = 1:n
        typ = methods(i);
        disp(typ)
        outputImage = getImage(filtKSMat{i,1}, filtKSMat{i,2}, filtKSMat{i,3}); % rename this to the varname
        outputImage = outputImage(:,:,1);
        matname = strcat('slice2_image_', string(typ), '.mat'); % slice2_image_badData
        save(matname, 'outputImage');
    end
    


    %%%%%%%%%% TESTING FUSE/RESTORE or RESTORE/FUSE %%%%%%%%%
    % using the splices, compress together 
    % get the channels for baseline data
    RF = restore_fuse(filtKSMat{1,1}, filtKSMat{1,2}, filtKSMat{1,3});
    FR = fuse_restore(filtKSMat{1,1}, filtKSMat{1,2}, filtKSMat{1,3});
    
    % display to check working
    %% plotting scripts
    close all
    figure(1); 
    imagesc(RF);
    title("Restore/Fuse")
    axis image, 
    colormap gray;
    axis off
    
    figure(2); 
    imagesc(FR);
    title("Fuse/Restore")
    axis image, 
    colormap gray;
    axis off

    % get evaluation critera for RF and FR
    load('slice2_image_goodData.mat')
    goodData_Image = outputImage;
    psnr_RF = calcPSNR(goodData_Image, RF);
    ssim_RF = calcSSIM(goodData_Image, RF);
    nmse_RF = calcNMSE(goodData_Image, RF);
    psnr_FR = calcPSNR(goodData_Image, FR);
    ssim_FR = calcSSIM(goodData_Image, FR);
    nmse_FR = calcNMSE(goodData_Image, FR);
    
    %%%%%%%%%%%%%%%% EVALUATE %%%%%%%%%%%%%%%%
    % evaluate the filtered vs. good image here for all filters
    % store output in an array reporting the scores
    % for all filters
        % psnr = calcPSNR(good,filt) -- on image data
        % ssim = calcSSIM(good,filt) -- on image data
        % nmse = calcNMSE(good, filt) -- on image data
    % return all evaluation metrics in a summary

    % pull in the good matrix
    load('slice2_image_goodData.mat')
    goodData_Image = outputImage;
    
    % pull in the first filtered matrix
    load('slice2_image_badData.mat');
    badData_Image = outputImage;

    % create array of list of files to iterate through
    cur_matFiles = dir('slice2_image_*');
    class(cur_matFiles)
    evalResults = [0,0,0];
    % iterate through all the files, calculate metrics, populate matrix
    for i = 1:length(cur_matFiles)
        disp(cur_matFiles(i).name)
        load(string(cur_matFiles(i).name));
        filtered_Image = outputImage;
        psnr_val = calcPSNR(goodData_Image, filtered_Image);
        ssim_val = calcSSIM(goodData_Image, filtered_Image);
        nmse_val = calcNMSE(goodData_Image, filtered_Image);
        evalResults(i,1) = psnr_val;
        evalResults(i,2) = ssim_val;
        evalResults(i,3) = nmse_val;
    end
    
    % results follow the alphabetical order for slice2_image_ in this dir
    disp('   PSNR       SSIM      NMSE')
    disp(evalResults)

    save("slice2_eval.mat", 'evalResults');

    %%%%%%%%%%%%%%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%

    %%% CREATING IMAGE FROM K SPACE %%%
    function RF = restore_fuse(ch1, ch2, ch3)
        % restore
       Data_img(:,:,1) = ifftshift(ifft2(ch1),1);
       Data_img(:,:,2) = ifftshift(ifft2(ch2),1);
       Data_img(:,:,3) = ifftshift(ifft2(ch3),1);
       % disp(Data_img(1,1,1));
       clear_comp = linspace(10,0.1,size(Data_img,2)).^2; 
       clear_matrix = repmat(clear_comp,[size(Data_img,1) 1]);

       % fuse
       eye_raw  = sqrt( abs(squeeze(Data_img(:,:,1))).^2 + ...
           abs(squeeze(Data_img(:,:,2))).^2 + ...
           abs(squeeze(Data_img(:,:,3))).^2).* clear_matrix; 
        % disp(eye_raw(1,1,1));
        crop_x = [128 + 60 : 348 - 33]; % crop coordinates
        eye_raw = eye_raw(crop_x, :);
        eye_visualize = reshape(squeeze(eye_raw(:,:)),[128 128]); 
        std_within = 0.9925; 
        [aa, val] = hist(eye_visualize(:),linspace(0,max(...
                                    eye_visualize(:)),1000));
        thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
        RF = uint16(eye_visualize * 65536 / thresh);
    end

    %%% CREATING IMAGE FROM K SPACE %%%
    function FusRes = fuse_restore(ch1, ch2, ch3)

       clear_comp = linspace(10,0.1,size(ch1,2)).^2;
       clear_matrix = repmat(clear_comp,[size(ch1,1) 1]);
       
       % does not work because you lose the phase info
       %eye_fuse  = sqrt( abs(squeeze(ch1(:,:))).^2 + ...
           %abs(squeeze(ch2(:,:))).^2 + ...
           %abs(squeeze(ch3(:,:))).^2).* clear_matrix; 
       
       eye_fuse= ((ch1 + ch2 + ch3)./3).* clear_matrix;

       restore = ifftshift(ifft2(eye_fuse),1);
       restore = abs(restore);
     
       crop_x = [128 + 60 : 348 - 33]; % crop coordinates
       restore = restore(crop_x, :);
       eye_visualize = reshape(squeeze(restore(:,:)),[128 128]); 
       std_within = 0.9925; 
       [aa, val] = hist(eye_visualize(:),linspace(0,max(...
                                    eye_visualize(:)),1000));
       thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
       FusRes = uint16(eye_visualize * 65536 / thresh);
        
    end

    %%% CREATING IMAGE FROM K SPACE %%%
    function eye_visualize = getImage(ch1, ch2, ch3)
       Data_img(:,:,1) = ifftshift(ifft2(ch1),1);
       Data_img(:,:,2) = ifftshift(ifft2(ch2),1);
       Data_img(:,:,3) = ifftshift(ifft2(ch3),1);
       clear_comp = linspace(10,0.1,size(Data_img,2)).^2; 
       clear_matrix = repmat(clear_comp,[size(Data_img,1) 1]);
       eye_raw  = sqrt( abs(squeeze(Data_img(:,:,1))).^2 + ...
           abs(squeeze(Data_img(:,:,2))).^2 + ...
           abs(squeeze(Data_img(:,:,3))).^2).* clear_matrix; 
        crop_x = [128 + 60 : 348 - 33]; % crop coordinates
        eye_raw = eye_raw(crop_x, :);
        eye_visualize = reshape(squeeze(eye_raw(:,:)),[128 128]); 
        std_within = 0.9925; 
        [aa, val] = hist(eye_visualize(:),linspace(0,max(...
                                    eye_visualize(:)),1000));
        thresh = val(find(cumsum(aa)/sum(aa) > std_within,1,'first'));
        eye_visualize = uint16(eye_visualize * 65536 / thresh);

    end
    
    %%% CALC PSNR %%%
    % higher values are better, peak signal to noise ratio
    function psnrValue = calcPSNR(good, filt)
        psnrValue = psnr(filt,good);
    end

    %%% CALC SSIM %%%
    % accounts for structural info, luminance, contrast, texture simularity
    % close to 1 indicates good image
    function ssimValue = calcSSIM(good, filt)
        ssimValue = ssim(good, filt);
    end

    %%% CALC NMSE %%%
    % ideal to be as low as possible
    function nmseValue = calcNMSE(good, filt)
        mse = mean((good(:) - filt(:)).^2);
        normConst = mean(good(:).^2);
        nmseValue = mse / normConst;
    end

end
