function [T,x,y] = SL_detection(K, M)                          % input corrupted K-space image
    % M= size(K,2);
    K= abs(K);   
                                                            
    for i = 1:M-2                                           % compute the sum of the difference between adjacent columns 
       T(i)= sum(abs(K(:,i+1)-K(:,i) + K(:,i+1)-K(:,i+2))); 
    end

    T = [0, T];
    [x,y] = findpeaks(T,'MinPeakProminence',75);             % detect large differences between adjacent columns 
end
