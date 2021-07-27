% Illustration of svd applied in image processing
% @Coyright Phi Ha

clear all; close all; clc
%pkg load image

%reading and converting the image
inImage = imread('Baby_AJ.jpg');  % over memory in OCTAVE while using svd
% inImage = imread('CR72.jpg');   % works fine with OCTAVE but not 
% really useful

if size(size(inImage)) ~= 2
   inImageD = rgb2gray(inImage);
 else 
   inImageD = inImage;
end

inImageD = double(inImageD);
rank(inImageD)

% decomposing the image using SVD
[U,S,V]=svd(inImageD);

% Using different number of singular values (diagonal of S) to compress and
% reconstruct the image
dispEr = []; numSVals = [];

for N = 10:20:90
    % store the singular values in a temporary var
    C = S;
    % discard the diagonal values not required for compression
    C(N+1:end,:)=0;
    C(:,N+1:end)=0;

    % Construct an Image using the selected singular values
    D=U*C*V';

    % display and compute error
    figure;
    buffer = sprintf('Image output using %d singular values', N)
    imshow(uint8(D));
    title(buffer);
    error=sum(sum((inImageD-D).^2));

    % store vals for display
    dispEr = [dispEr; error];
    numSVals = [numSVals; N];
    
    if N == 80
      print -djpeg Baby_AJ1.jpg
    end
end

%% dislay the error graph
%figure; 
%title('Error in compression');
%plot(numSVals, dispEr);
%grid on
%xlabel('Number of Singular Values used');
%ylabel('Error between compress and original image');

