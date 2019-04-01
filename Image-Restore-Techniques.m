clear all;
close all;
clc;

% 1) Wiener Filter (Least Mean Square Filter)

I = imread('windmill.bmp');
figure(1)
imshow(I)
title('Image with no blur or noise');

h = fspecial('gaussian',[5 5], 5);  %δημιουργία φίλτρου gaussian [a b] = mask
Ifiltered = imfilter(I,h);  
figure(2)
imshow(Ifiltered)
title('Blurred image');

mean = 0;
variance = 0.001;
J = imnoise(Ifiltered,'gaussian',mean,variance); % WGN
figure(3)
imshow(J)
title('Blurred and noise image');

[x, y] = size(I);
H = fft2(h,x,y);
Hc = conj(H);
G = fft2(J);
Sn = variance;
Sf = abs(fftshift(fft2(J))).^2;

gamma = 1;
F = (Hc./(abs(H).^2 + gamma*100000000*(Sn./Sf))).*G;  
f = ifft2(F);

figure(4);
imshow(uint8(f));

mse1 = sqrt(sum(sum((double(I) - double(f)) .^ 2))) / (x * y)


% 2) Constrained Least Squares Restoration

I2 = imread('windmill.bmp');
figure(5)
imshow(I2)

h2 = fspecial('gaussian',[5 5], 5); 
Ifiltered2 = imfilter(I2,h2);  
figure(6)
imshow(Ifiltered2)

mean2 = 0;
variance2 = 0.005;
J2 = imnoise(Ifiltered2,'gaussian',mean2,variance2); % WGN
figure(7)
imshow(J2)

p = [0 -1 0;-1 4 -1;0 -1 0]; 

[A, B] = size(I2);
M = A+3-1;
N = B+3-1;
Pex = fft2(p,A,B);

G = fft2(J2,A,B);
g = ifft2(G);

H = fft2(h2,A,B);
Hstar = conj(H);

gamma = 1;
square_mean = mean2.^2;
n = norm((A - 1)*(B - 1)*(square_mean + variance2));
trm = 0;
a = (80*n^2)/100;

while(trm ~= 1)
    F = (Hstar./(abs(fftshift(H)).^2 + gamma*(abs(fftshift(Pex)).^2))).*G;
    f = ifft2(F);    
    r = G - H.*F;
    phi_gamma = norm(r).^2;
       
    if(phi_gamma > (n - a) & phi_gamma < (n + a))
        trm = 1;
    elseif(phi_gamma < (n - a))
        gamma = gamma + 0.1*gamma;
    elseif(phi_gamma > (n - a))
        gamma = gamma - 0.1*gamma;
    end;
   % gamma
end

figure(8);
imshow(uint8(f));

mse2 = sqrt(sum(sum((double(I2) - double(f)) .^ 2))) / (x * y)


