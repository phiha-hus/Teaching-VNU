clear all; close all; clc

%pkg load image
%pkg load imagesc

%iptsetpref('ImshowBorder','tight') 
L=imread('CR72.jpg');

imshow(L)

L1=L(:,:,1); 
L2=L(:,:,2);
L3=L(:,:,3);

I1=im2double(L1); I2=im2double(L2); I3=im2double(L3);

figure
imshow(I1)
figure
imshow(I2)
figure
imshow(I3)
[u1,s1,v1]=svd(I1);
[u2,s2,v2]=svd(I2);
[u3,s3,v3]=svd(I3);
C1=zeros(size(I1));
C2=zeros(size(I2));
C3=zeros(size(I3));
k=15;
for j=1:k
C1=C1+s1(j,j)*u1(:,j)*v1(:,j).';
end
for j=1:k
C2=C2+s2(j,j)*u2(:,j)*v2(:,j).';
end
for j=1:k
C3=C3+s3(j,j)*u3(:,j)*v3(:,j).';
end
figure
imshow(C1)

figure
imshow(C2)
figure
imshow(C3)
C1(k)=1;
C2(k)=1;
C3(k)=1;
R1=im2uint8(C1);
R2=im2uint8(C2);
R3=im2uint8(C3);
Q(:,:,1)=R1;
Q(:,:,2)=R2;
Q(:,:,3)=R3;
figure
imshow(Q,[])
