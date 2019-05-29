%% synthetic aperture for 64 scanlines the right one
clc;
clear all
close all
cd 'C:\Users\Valyria\Documents\Ultrasound'
load('anecoicCystData.mat')
v=1540;
no_lines = 128;
es=veraStrct.elementSpacingMM;
elementposmm=(-127/2*es:es:127/2*es)*10^-3;
focal_point_es = linspace(0,90*1e-3,2353);
num_elements = 128;
f0=veraStrct.frequencyMHz*1e6;%center frequency
fs1 = veraStrct.samplingRateMHz*1e6;
c = 1540;
data_pt = veraStrct.data(80:end,:,:);
depth = data_pt(:,1:128,65);%choosing lateral location at 65 - A-line 65
pixelintomm = (0.5*c)/fs1;
% min_depth = min(min(depth))*pixelintomm;
%depth_x = (0:size(depth,1)+min_depth);
pixnum = 1:1:length(depth);
pixmm = pixnum *pixelintomm;
elementposmm=(-127/2*es:es:127/2*es)*1e-3;
scanlines = 64; %for now
f0 = 0.0015;%transmit focus in m
p_depth_focus = 0.04;
% you use different elements for recieveing and transmitting
% for i = 4:1:124
%     p0(i) = (elementposmm(i-3)+elementposmm(i-2)+elementposmm(i-1)+elementposmm(i)+elementposmm(i+1)+elementposmm(i+2)+elementposmm(i+3)+elementposmm(i+4))/8;
% end
% for j = 1:3
%     p1(j) = (elementposmm(j)+elementposmm(j+1)+elementposmm(j+2)+elementposmm(j+3)+elementposmm(j+4))/4;
% end 
% for k = 125:128
%     p2(k) = (elementposmm(k-3)+elementposmm(k-2)+elementposmm(k-1)+elementposmm(k))/4;
% end
% pvals = [p1';p0(4:end)';p2(125:end)'];
% fprintf('the signal error iir is:%7.2f  \n',signal_error_iir)
zf = f0;
dp = linspace(0,0.09,2353);

for i = 1:1:128
    for j = 1:1:2353
        time_array(i,j) = j/fs1;
    end
end
p0=(v*time_array(1,:))/2;
p0 = p0';
%ti = find(dp==p_depth_focus );
ts = (1:1:2353)/fs1;
%time = 2*t(1047);%set it manually to see where the pdepth is this is in s - roundtrip time so multiply by 2
v = 1540;
xr = es*1e-3;%in m
xf = xr;%in m
%xr = elementposmm;
xfp3 = 0;
ele256beam = ((-256/2)*es:es:(256/2*es))*1e-3;
%elementposmm = 1*es:es:128*es;
%for k = 1%elements that would receive the beam

   %for j = 1:128%number of elements that would transmit the beam
   %j = 0;
%for m = 1:length(ele256beam)
%        cd 'C:\Users\Valyria\Documents\Ultrasound'
%        p0 = (v*time_array(1,:))/2+(m*xf);
for m = 0:127
    for j = 1:length(elementposmm)
     for i =1:1:2353 %yposiyion for p00
      for k = 1:1:128% xposition for p0
       
       
       %for k = 1:1:128
            if p0(i) < p_depth_focus
              % image_scanline(i,k)= ((p_depth_focus-sqrt(((dp(i)-p_depth_focus).^2)+elementposmm(j)^2))/v+(sqrt(dp(i)^2)+(elementposmm(k)-elementposmm(j)).^2)/v);
                image_scanline(i,j,k)= ((p_depth_focus-sqrt(((p0(i)-p_depth_focus).^2)+(0)^2))/v+(sqrt(p0(i)^2+((ele256beam(j))-((ele256beam(k+m)))).^2)/v));
%                 image_scanline_t(i,k,j) = (p_depth_focus-sqrt(((p0(i)-p_depth_focus).^2)+0^2))/v;
%                 ir(i,k) = (sqrt(p0(i)^2+(elementposmm(k)-elementposmm(65)).^2)/v);
            elseif p0(i) > p_depth_focus
                 image_scanline(i,j,k)= ((p_depth_focus+sqrt(((p0(i)-p_depth_focus).^2)+(0)^2))/v+(sqrt(p0(i)^2+((ele256beam(j))-((ele256beam(k+m)))).^2)/v));
%                  image_scanline_t(i,k) = (p_depth_focus+sqrt(((p0(i)-p_depth_focus).^2)+0^2))/v;
%                 ir(i,k) = (sqrt(p0(i)^2+(elementposmm(k)-elementposmm(65)).^2)/v);
               %image_scanline(i,k)= ((p_depth_focus+sqrt(((dp(i)-p_depth_focus).^2)+elementposmm(j)^2))/v+(sqrt(dp(i)^2)+(elementposmm(k)-elementposmm(j)).^2)/v);
               %0 is elementposmm(m) later 
    
            %image_scanline(i,j,k)= ((2*dp(i))/v)-((p_depth_focus-sqrt(((dp(i)-p_depth_focus).^2)+0^2))/v+(sqrt(dp(i)^2)+(elementposmm(j)-elementposmm(k)).^2)/v);
%             else 
%                 image_scanline(i,j,k)= (2*dp(i))/v-((p_depth_focus+sqrt(((dp(i)-p_depth_focus).^2)+0^2))/v+(sqrt(dp(i)^2)+(elementposmm(j)-elementposmm(k)).^2)/v);
%         s = 'scanline';
%         num = num2str(j);
        %tt = 2*(p_depth_focus/v);
        %filename = strcat(s,num);
        %image_scanline3D(i,j) =((dp(i)-sqrt(dp(i)-(v*tt/2).^2+xr^2))/v+(sqrt((v*tt/2).^2+(elementposmm(j)-elementposmm(65)).^2))/v);
%         save(filename,'image_scanline3D');
%         disp(i)
        %image_scanline(i,j) = tt-((dp(i)-sqrt(dp(i)-(v*tt/2).^2+xr^2))/v+(sqrt((v*tt/2).^2+(elementposmm(65)-elementposmm(j)).^2))/v);
        %((v*tt/2).^2 +(elementposmm(i)-elementposmm(j)).^2)/v));
            %end      %end
    end
   end
  end
end
   s = 'scanline';
   num = num2str(m+1);
   filename = strcat(s,num);
   cd 'C:\Users\Valyria\Documents\Ultrasound\scanlines'
   save(filename,'image_scanline')
end

%ImageData = zeros(2353,128,128);
% delay3 = squeeze(image_scanline(:,:,3));
% delay2 = squeeze(image_scanline(:,:,2));
disp ('done')
%%
%load('scanline0.mat')
% timeipix =image_scanline *(fs1);
clc
for i = 1:1:128
    for j = 1:1:2353
        time_array(i,j) = j/fs1;
    end
end
x = time_array';
pl = 128;
%k = 2;
for k = 3:130
    cd 'C:\Users\Valyria\Documents\Ultrasound\scanlines'
    Files=dir('*.*');
    FileNames=Files(k).name;
    load(FileNames)  
    
    for i = 1:pl %Aline
        DataforALine = squeeze(data_pt(:,:,i));
        imagesca = squeeze(image_scanline(:,:,i));
        
            for j = 1:128%element
              
              ChDataDelay(:,j) = interp1(x(:,j),DataforALine(:,j),((imagesca(:,j))));
              ImageData(:,j,i) = ChDataDelay(:,j);
            end
    end
    s = 'interpdata';
   num = num2str(k-2);
   filename = strcat(s,num);
   cd 'C:\Users\Valyria\Documents\Ultrasound\interpdata'
   save(filename,'ImageData')
end
disp('done')



%% right way?!
clc
for k = 3:130
        cd 'C:\Users\Valyria\Documents\Ultrasound\interpdata'
        Files=dir('*.*');
        FileNames=Files(k).name;
        load(FileNames) 
        sumim = sum(ImageData,2);
        sumsq = squeeze(sumim);
        s = 'summed_data';
        num = num2str(k-2);
        filename = strcat(s,num);
        cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
        save(filename,'sumsq')
end
%% 
fin_image = zeros(2353,256);
cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
load('summed_data1.mat')
fin_image(:,1:128) = sumsq(:,:);
for k = 4:130
        cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
        Files=dir('*.*');
        FileNames=Files(k).name;
        load(FileNames) 
        fin_image(:,k-2:k-2+127) = fin_image(:,k-2:k-2+127)+sumsq(:,:);
        
end
sumH = (20*log10(abs(hilbert(fin_image(2:end,1:end)')))); 
sumH = -max(max(sumH))+sumH;
imagesc(sumH',[-50,0]);
%axis image
colormap(gray)
%%

% im_data = zeros(2353,128,256);
% for i = 2:128
%             im_data(:,i,i:i+127)=im_data(:,i-1,i:i+127)+ImageData(:,i,1:128);
% end
%    sumIm = sum(im_data,2);
%    sumIm2 = squeeze(sumIm);
%    im_scan =(20*log10(abs(hilbert(sumIm2')))); 
%    im_scan = -max(max(im_scan))+im_scan;
% 
% 
%     imagesc(elementposmm,pixmm,im_scan',[-60,0]);
%     axis image
%     colormap(gray)
    %% just for test purposes
    
    image1 = sumsq;
    image2 = zeros(2353,129);
    image3 = zeros(2353,130);
    image4 = zeros(2353,131);
    image5 = zeros(2353,132);
    image6 = zeros(2353,133);
    image7 = zeros(2353,134);
    image8 = zeros(2353,135);
    image9 = zeros(2353,136);
    image10 = zeros(2353,137);
    
    imagef2 = zeros(2353,137);
%     imagefinal2(:,1:128) = image1;
%     image2(:,1:128) = image1;
    %imagefinal = image3(:,3:130)+image2(:,2:129)+image1(:,:)+image4(:,4:131)+image5(:,5:132)+image6(:,6:133)+image7(:,7:134)+image8(:,8:135)+image9(:,9:136)+image10(:,10:137);
    %imagef2(:,:) =imagefinal2(:,2:129)+image1(:,:);% image2(:,:)+image3(:,3:130);
    imagef2(:,1:128) = image1(:,:);
for i = 2:10 
     imagef2(:,i:i+127) = imagef2(:,i:i+127)+image1(:,1:128);
end
    
    



    subplot (121)
    im_scan2 =(20*log10(abs(hilbert(imagef2')))); 
    im_scan2 = -max(max(im_scan2))+im_scan2;
    imagesc(elementposmm,pixmm,im_scan2',[-50,0]);
    axis image
    colormap(gray)
    title('3 sum imaged')
    
    subplot (122)
    sumim = sum(ImageData,2);
    sumsq = squeeze(sumim);
    sumH = (20*log10(abs(hilbert(sumsq')))); 
    sumH = -max(max(sumH))+sumH;
    imagesc(elementposmm,pixmm,sumH',[-50,0]);
    axis image
    colormap(gray)
    title('single image')
    %%
   
%         sumNewImage = sum(ImageData,2);
%         sumNewImage = squeeze(sumNewImage);
    
%         s = 'sumImage';
%         num = num2str(k-2);
%         filename = strcat(s,num);
%         cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
%         save(filename,'sumNewImage')
%end
disp('done')
%% wrong way
clc
element_256posmm = (-256/2*es:es:256/2*es)*1e-3;
total_scan = zeros(2353,256);
total_scan(:,1:128) = sumNewImage(:,1:128);

for i = 2:128
     %for j = 1:128
        total_scan(:,i:i+127) = total_scan(:,i:i+127)+sumNewImage(:,:);
     %end 
 end
% cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
% load('sumImage1.mat')
%(:,1:128)=sumNewImage(:,:);
%total_scan(:,1:128) = ImageData(:,:,1);
% for k = 4:130
%         cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
%         Files=dir('*.*');
%         FileNames=Files(k).name;
%  for i = 1:128 
% for i = 1:128
%     for j = 1:128
%     
%             total_scan(:,j:j+127,i)=total_scan(:,j:j+127,i)+ImageData(:,:,i);
%         
%     end
% end
% totalscan = totalscan(:,i:i+128)+totalscan(:,:,i);
% xm = ImageData(:,:,92);


totaltime = total_scan;
im_scan =(20*log10(abs(hilbert(totaltime')))); 
im_scan = -max(max(im_scan))+im_scan;


imagesc(elementposmm,pixmm,im_scan',[-50,0]);
axis image
colormap(gray)
%  end
%  
 
disp('done')
%% wrong mask


fnum = 2;
depth = linspace(1,90,2353);
d = (depth/fnum)*10^-3;
%elementposmm=(1*es:es:128*es);
elementposmm=(-230/2*es:es:230/2*es)*10^-3;
for j = 1:1:2353
    
    
        elementp(j) = round((elementposmm(125)+(d(j)*(1/3)))*(10^3/es)+124);
        elementn(j) = round((elementposmm(125)-(d(j)*(2/3)))*(10^3/es)+124);
    
    
end
%d in is at 1308

dines = round(d(1:end)*(10^3/es));
dines2 = find(dines==230);


% elementnpos = round(((elementn*10^3)/es));
% elementppos = round(((elementp*10^3)/es));
bwmask = ones(2353,256);
bwmask(129,1) =0 ;
for m = 1:1:230
    for n = 1:1:dines2(end)
            valp = n;
            elementinn = 115-(dines(n)/2);
            elementinp = 115+(dines(n)/2);
%             disp(elementinp)
%             disp(elementinn)
            if elementinn == 0
                elementinn = 1;
            end
            bwmask(valp,1:round(elementinn)) = 0;
             bwmask(valp,round(elementinp):end) = 0;
            
%             valp = elementp(n);
%             valn = elementn(n);
%             valpn = elementn<=valn;
%             valpp =  elementp>=valp;
%             xm = find(valpn==1);
%             xn = find(valpp==1);
%             bwmask(xm,valn) = 1;
%             bwmask(xn,valp) = 1;
%     bwmask(elementppos>=valp,i) = 0;
    end 
end

% for i = 1:1:1311
%     mn = 
%     bmask(i,elementinn(i)<elementinn(i)) = 1;
%     bmask(i,elementinp(i)>elementinp(i)) = 1;
% end
    
% bwmask(valp,round(elementinn)) = 1;
% bwmask(valp,round(elementinp)) = 1;
figure;
imagesc(bwmask)
colormap gray
%% right mask


% for i = 1:1:2353
%     vals = elementnpos< round(elementnpos(i));
%     bwmask(i,vals) = 0;
% end



%subplot(121)

%elementposmm=(-127/2*es:es:127/2*es)*10^-3;
% masksum = sum(masking,2);
% mask_sq = squeeze(masksum);

bm = ones(2535,256);
mn = im_scan';
first_half = mn(1:2353,1:100)<=-39;
np = zeros(2353,256);
np = (first_half==1);
for i = 1:2353
    for j = 1:100
        if np(i,j) == 1
            bm(i,j) = 0;
        else
            bm(i,j) = 1;
        end
    end
end
second_half = mn(1:2353,100:256)<=-39;
hj = (second_half==1);
gg = size(hj);
gk = gg(2);
for i = 1:2353
    for j = 1:gk
        if hj(i,j) == 1
            bm(i,100+j) = 0;
        else
            bm(i,100+j) = 1;
        end
    end
end


% imagesc(bm)
% colormap gray
for j = 1:1:2353
    for i = 1:1:256
            masking(j,i) = (sumIm2(j,i).*bm(j,i));
            disp(j)
            
    end
end
mask_sq2 = (20*log10(abs(hilbert(masking(2:end,:)'))));
bmode_mask = -max(max(mask_sq2))+mask_sq2;
imagesc(element_256posmm,pixmm,bmode_mask',[-60,0])
colormap gray
axis image
    

% tt = totaltime.*bwmask(:,:);
% im_scan2 =(20*log10(abs(hilbert(tt')))); 
% im_scan2 = -max(max(im_scan2))+im_scan2;
% 
% 
% imagesc(elementposmm,pixmm,im_scan2',[-70,0]);
% %axis image
% colormap(gray)
            


