function R = Autocorrelation(X,dimx,dimy,p)     
    %patch size p can only be a common divisor of 512 and 128
    x= zeros(dimx/p,1)+p;                       % define the x dimensions of the k-subspaces
    y= zeros(dimy/p,1)+p;                       % define the y dimensions of the k-subspaces
    Y= mat2cell(X,x,y);                         % create the square submatrices
    
    for m = 1:length(x)
        for n= 1:length(y)
            T = sum(Rh(cell2mat(Y(m,n)))+ Rh(transpose(cell2mat(Y(m,n)))),1); % compute the sum of total vertical and horizontal correlations
            V = flip(1:p)*2*p;                                                % define the proper averaging coëfficients 
            R{m,n} = bsxfun(@rdivide,T,V);                                    % divide the total sum of correlations by the correct averaging constant  
    
        end
    end
    %%%%%%% GET HORIZONAL AUTOCORR %%%%%%%
    function Rh= Rh(A)  % determines horizontal autocorrelations of a square matrix                      
        p = size(A,1);                               
            for i = 1:p                             
                    for k= 0:p-1
                        for j = 1:p
                             if j+k>p 
                                break
                             else
                                a(j)= A(i,j)*conj(A(i,j+k));
                         K(i,k+1)=sum(a);
                             end 
                        end
                     a = [];
                    end  
            end
        Rh=conj(K);
     end
end           
    
                   
