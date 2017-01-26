clear;
close all;
clc;

%%read and store multiple image frames
file_path = './EnterExitCrossingPaths2cor/';
img_path_list = dir(strcat(file_path,'*.jpg'));
img_num = length(img_path_list); %% the number of image frames
img = cell(1,img_num); %% used to store original image frames
img_gray = cell(1,img_num); %%used to store gray values of image frames
for i=1:img_num
    image_name=img_path_list(i).name;
    image_name=strcat(file_path,image_name);
    img{i}=imread(image_name);
    img_gray{i}=rgb2gray(img{i});%% make image frames grayscale
end

h3=fspecial('average',[3,3]);
h5=fspecial('average',[5,5]);
G1=fspecial('Gaussian',[5,5],1.0);
G2=fspecial('Gaussian',[7,7],1.4);
G3=fspecial('Gaussian',[9,9],1.8);
I_h3=cell(1,img_num);
I_h5=cell(1,img_num);
I_G1=cell(1,img_num);
I_G2=cell(1,img_num);
I_G3=cell(1,img_num);
for i=1:img_num
    %%spatial smoothing filter with 3X3 box filter
    I_h3{i}=imfilter(img_gray{i},h3,'replicate');
    %%spatial smoothing filter with 5X5 box filter
    I_h5{i}=imfilter(img_gray{i},h5,'replicate');
    %%spatial smoothing filter with 2D Gaussian filters with standard deviation ssigma=1.0
    I_G1{i}=imfilter(img_gray{i},G1,'replicate');
    %%spatial smoothing filter with 2D Gaussian filters with standard deviation ssigma=1.4
    I_G2{i}=imfilter(img_gray{i},G2,'replicate');
    %%spatial smoothing filter with 2D Gaussian filters with standard deviation ssigma=1.8
    I_G3{i}=imfilter(img_gray{i},G3,'replicate');
end

len=size(img{1},2); %%the length value of the matrix correspond to one image frame
wid=size(img{1},1); %%the width value of the matrix correspond to one image frame
sf_G1=cell(wid,len);%%sf_G1 is a?wid x len?matrix, each entry is also a (1 x img_num) matrix
%%store each corresponding gray values of multiple image frames filtered by a 2D Gaussian filter
for i=1:wid
    for j=1:len
        for k=1:img_num
            sf_G1{i,j}(k)=I_G1{k}(i,j);
        end
    end
end
gx1_G1=sf_G1;
gx2_G1=sf_G1;
gx3_G1=sf_G1;
sim_f=[-0.5,0,0.5]; %%create a simple 0.5[-1, 0, 1] filter
t1=1.0; %%value of sigma t1
t2=1.4; %%value of sigma t2
t3=2.6; %%value of sigma t3
gx1_mask=zeros(1,5);%range 5*t1
gx2_mask=zeros(1,7);%range 5*t2
gx3_mask=zeros(1,13);%range 5*t3

%%a 1D derivative of Gaussian filter with standard devition of 1.0
for x=-2:2
    gx1_mask(x+3)=-x/(t1*t1)*exp(-(x*x)/(2*t1*t1));
end
%%a 1D derivative of Gaussian filter with standard devition of 1.4
for x=-3:3
    gx2_mask(x+4)=-x/(t2*t2)*exp(-(x*x)/(2*t2*t2));
end
%%a 1D derivative of Gaussian filter with standard devition of 1.8
for x=-6:6
    gx3_mask(x+7)=-x/(t3*t3)*exp(-(x*x)/(2*t3*t3));
end

for i=1:wid
    for j=1:len
        sf_G1{i,j}=imfilter(double(sf_G1{i,j}),sim_f,'replicate');
        gx1_G1{i,j}=imfilter(double(gx1_G1{i,j}),gx1_mask,'replicate');
        gx2_G1{i,j}=imfilter(double(gx2_G1{i,j}),gx2_mask,'replicate');
        gx3_G1{i,j}=imfilter(double(gx3_G1{i,j}),gx3_mask,'replicate');
    end
end

thr_sf_G1=cell(1,img_num);
thr_gx1_G1=cell(1,img_num);
thr_gx2_G1=cell(1,img_num);
thr_gx3_G1=cell(1,img_num);
%store gray values of each image frames filtered by different temporal derivative filter
for k=1:img_num
    for i=1:wid
        for j=1:len 
            thr_sf_G1{k}(i,j)=sf_G1{i,j}(k);
            thr_gx1_G1{k}(i,j)=gx1_G1{i,j}(k);
            thr_gx2_G1{k}(i,j)=gx2_G1{i,j}(k);
            thr_gx3_G1{k}(i,j)=gx3_G1{i,j}(k);
        end
    end
end

%calculate the threshold
thres1=cell(1,img_num);
thres2=cell(1,img_num);
thres3=cell(1,img_num);
thres4=cell(1,img_num);
for i=1:img_num
    thres1{i}=(std(double(thr_sf_G1{i}(:)))+7)*ones(wid,len);
    thres2{i}=(std(double(thr_gx1_G1{i}(:)))+7)*ones(wid,len);
    thres3{i}=(std(double(thr_gx2_G1{i}(:)))+7)*ones(wid,len);
    thres4{i}=(std(double(thr_gx3_G1{i}(:)))+7)*ones(wid,len);
end

final1_img=cell(1,img_num);%%store the final mask after thresholding
final2_img=cell(1,img_num);
final3_img=cell(1,img_num);
final4_img=cell(1,img_num);
%Convert image to binary image by thresholding
for i=1:img_num
    thr_sf_G1{i}=im2bw(abs(thr_sf_G1{i})-thres1{i},1/255);
    thr_gx1_G1{i}=im2bw(abs(thr_gx1_G1{i})-thres2{i},1/255);
    thr_gx2_G1{i}=im2bw(abs(thr_gx2_G1{i})-thres3{i},1/255);
    thr_gx3_G1{i}=im2bw(abs(thr_gx3_G1{i})-thres4{i},1/255);
    %Combine the mask with the original frame
    final1_img{i}=double(img_gray{i}).*double(thr_sf_G1{i});
    final2_img{i}=double(img_gray{i}).*double(thr_gx1_G1{i});
    final3_img{i}=double(img_gray{i}).*double(thr_gx2_G1{i});
    final4_img{i}=double(img_gray{i}).*double(thr_gx3_G1{i});
end

%display the result of images filtered by a simple[-0.5 0 0.5] filter
for i=1:img_num
    imshow(final1_img{i});
    pause(0.05);
end
%display the result of images filtered by a 1D derivative 
%of Gaussian filter with standard devition of 1.0
for i=1:img_num
    imshow(final2_img{i});
    pause(0.1);
end
%display the result of images filtered by a 1D derivative
%of Gaussian filter with standard devition of 1.4
for i=1:img_num
    imshow(final3_img{i});
    pause(0.05);
end
%display the result of images filtered by a 1D derivative 
%of Gaussian filter with standard devition of 1.8
for i=1:img_num
    imshow(final4_img{i});
    pause(0.05);
end


 



