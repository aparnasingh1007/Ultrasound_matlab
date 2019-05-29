%% synthetic aperture for 64 scanlines the right one
clc;
clear all
close all
cd 'C:\Users\Valyria\Documents\Ultrasound'
load('pointTargetData.mat')
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
elementposmm=(-128/2*es:es:128/2*es)*1e-3;
elementposmm = elementposmm(1:128);
elementposmm(64) = 0;
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
ele256beam = (-256/2)*es:es:256/2*es;
%elementposmm = 1*es:es:128*es;
%for k = 1%elements that would receive the beam

   %for j = 1:128%number of elements that would transmit the beam
   %j = 0;
%for m = 1:length(ele256beam)
%        cd 'C:\Users\Valyria\Documents\Ultrasound'
%        p0 = (v*time_array(1,:))/2+(m*xf);
    for j = 1:length(elementposmm)
     for i =1:1:2353 %yposiyion for p00
      for k = 1:length(elementposmm)% xposition for p0
       
       
       %for k = 1:1:128
            if p0(i) < p_depth_focus
              % image_scanline(i,k)= ((p_depth_focus-sqrt(((dp(i)-p_depth_focus).^2)+elementposmm(j)^2))/v+(sqrt(dp(i)^2)+(elementposmm(k)-elementposmm(j)).^2)/v);
                image_scanline(i,j,k)= ((p_depth_focus-sqrt(((p0(i)-p_depth_focus).^2)+(es*1e-3)^2))/v)+((sqrt(p0(i)^2+((elementposmm(k))-((elementposmm(j)))).^2))/v);
%                 image_scanline_t(i,k,j) = (p_depth_focus-sqrt(((p0(i)-p_depth_focus).^2)+0^2))/v;
%                 ir(i,k) = (sqrt(p0(i)^2+(elementposmm(k)-elementposmm(65)).^2)/v);
            elseif p0(i) > p_depth_focus
                 image_scanline(i,j,k)= ((p_depth_focus+sqrt(((p0(i)-p_depth_focus).^2)+(es*1e-3)^2))/v)+((sqrt(p0(i)^2+((elementposmm(k))-((elementposmm(j)))).^2))/v);
%                
    end
   end
  end
end
%    s = 'scanline';
%    num = num2str(m);
%    filename = strcat(s,num);
%    cd 'C:\Users\Valyria\Documents\Ultrasound\scanlines'
%    save(filename,'image_scanline')
%  

%ImageData = zeros(2353,128,128);
% delay3 = squeeze(image_scanline(:,:,3));
% delay2 = squeeze(image_scanline(:,:,2));
disp ('done')
%%
%load('scanline0.mat')
% timeipix =image_scanline *(fs1);
clc
clear time_array
for k = 1:128
    for i = 1:1:128
        for j = 1:1:2353
         time_array(j,i,k) = j/fs1;
        end
    end
end
x = time_array;
% %%
% for i = 1:1:128
%     for j = 1:1:2353
%         time_array(i,j) = j/fs1;
%     end
% end
% x = time_array';
pl = 128;
%k = 2;
%for k = 3:130
%     cd 'C:\Users\Valyria\Documents\Ultrasound\scanlines'
%     Files=dir('*.*');
%     FileNames=Files(k).name;
%     load(FileNames)  
%i in imagesc is receive, i in dataforaline is transmit
    
for k = 1:128  
    %DataforALine = squeeze(data_pt(:,:,k));
    for i = 1:pl %Aline
      % imagesca = squeeze(image_scanline(:,:,i));
        
            for j = 1:128%element
              
              ChDataDelay(:,j,i) = interp1(x(:,j,i),data_pt(:,j,k),((image_scanline(:,j,i))),'linear');
              %ChDataDelay(:,j) = interp1(x(:,j),DataforALine(:,j),((imagesca(:,j))),'linear');
              ImageData(:,j,i) = ChDataDelay(:,j,i);
            end
    end
    %store that in a variable indexed by the outermost index loop
    
    s = 'interpdata';
   num = num2str(k);
   filename = strcat(s,num);
   cd 'C:\Users\Valyria\Documents\Ultrasound\interpdata'
   save(filename,'ImageData')
   disp(k)
end
 disp('done')
% sumdat = zeros(2353,256);
% %load('interpdata1.mat')
% summ = squeeze(sum(ImageData,2));
% sumdat(:,1:128) = summ(:,:);
% sumH = (20*log10(abs(hilbert(sumdat(2:2308,1:128))))); 
% sumH = -max(max(sumH))+sumH;

% subplot(121)
% imagesc(elementposmm,pixmm,sumH,[-70,0]);
% % axis image
%  colormap(gray)


%% right way?!
clc
% sumdat = zeros(2353,256);
% load('interpdata1.mat')
% summ = squeeze(sum(ImageData,2));
% sumdat(:,1:128) = summ(:,:);
for k = 1:128
        cd 'C:\Users\Valyria\Documents\Ultrasound\interpdata'
        %Files=dir('*.*');
        s = 'interpdata';
        num = num2str(k);
        FileNames = strcat(s,num);
%         Files=dir('*.*');
%         FileNames=Files(i+1).name;
        load(FileNames) 
        
        sumim = sum(ImageData,2);
        sumsq = squeeze(sumim);
        %sumdat(:,k-2:k-2+127) = sumdat(:,k-2:k-2+127)+sumsq(:,:);
        s = 'sumdata';
        num = num2str(k);
        filename = strcat(s,num);
        cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
        save(filename,'sumsq')
        disp(k)
end
%%
clc

sumdat = zeros(2353,256);
% load('sumdata1.mat')
%summ = squeeze(sum(ImageData,2));
%sumdat(:,1:128) = sumsq(:,:);
% sumH = (20*log10(abs(hilbert(sumdat(2:2308,1:242))))); 
% sumH = -max(max(sumH))+sumH;
% 
% %subplot(121)
% imagesc(elementposmm,pixmm,sumH,[-50,0]);
% axis image
% colormap(gray)
% title('one beam')
for i = 1:128
        cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
        s = 'sumdata';
        num = num2str(i);
        FileNames = strcat(s,num);
%         Files=dir('*.*');
%         FileNames=Files(i+1).name;
        load(FileNames) 
        sumdat(:,i:i+127) = sumdat(:,i:i+127)+sumsq(:,:);
         
end
%    sumIm = sum(im_data,2);
%    sumIm2 = squeeze(sumIm);
ele_256posmm = (-255/2*es:es:255/2*es)*1e-3;
   im_scan =(20*log10(abs(hilbert(sumdat(2:2308,:))))); 
   im_scan = -max(max(im_scan))+im_scan;
% 
%     subplot(122)
    figure;
    subplot(121)
    imagesc(ele_256posmm,pixmm,im_scan,[-50,0]);
    axis image
    colormap(gray)
    title('Unmasked data')
%     title('all 128 beams')
%    
%imagesc(im_scan,[-50,0]);colormap gray
  







    
    
%% right mask
close all

% sumdat = zeros(2353,256);
cd 'C:\Users\Valyria\Documents\Ultrasound\interpdata'
load('interpdata1.mat')
summ = squeeze(sum(ImageData,2));
sumdat(:,1:128) = summ(:,:);
sumH = (20*log10(abs(hilbert(sumdat(2:2308,1:128))))); 
sumH = -max(max(sumH))+sumH;
figure;
imagesc(sumH,[-70,0]);
%axis image
colormap(gray)


fnum = 2;
depth = linspace(1,90,2353);
d = (depth/fnum)*10^-3;
dines = round(d(1:end)*(10^3/es));
dines2 = find(dines==53);

%subplot(121)

bm = ones(2353,128);
%mn = im_scan;
%bm(200,1:30) = 0;

bm(200:250,1:26) =0 ;
bm(251:300,1:31) = 0;
bm(301:350,1:34) = 0;
bm(351:400,1:36) = 0;
bm(401:451,1:40) = 0;
bm(451:500,1:41) = 0;
bm(501:550,1:43) = 0;
bm(551:600,1:45) = 0;
bm(601:650,1:47) = 0;
bm(651:700,1:49) = 0;
bm(700:751,1:51) = 0;
bm(751:800,1:53) =0 ;
bm(800:850,1:55) = 0;
bm(851:900,1:57) = 0;
bm(901:950,1:57) = 0;
bm(951:1000,1:57) = 0;
bm(1001:1050,1:57) = 0;
bm(1051:1100,1:55) = 0;
bm(1100:1150,1:54) = 0;
bm(1151:1200,1:52) = 0;
bm(1201:1250,1:50) = 0;
bm(1251:1300,1:48) = 0;
bm(1301:1350,1:46) = 0;
bm(1351:1400,1:44) = 0;
bm(1401:1450,1:42) = 0;
bm(1451:1500,1:40) = 0;
bm(1501:1550,1:38) = 0;
bm(1551:1600,1:36) = 0;
bm(1601:1650,1:34) = 0;
bm(1651:1700,1:32) = 0;
bm(1701:1750,1:30) = 0;
bm(1751:1800,1:28) = 0;
bm(1801:1900,1:26) = 0;
bm(1901:1950,1:24) = 0;
bm(1951:2000,1:22) = 0;
bm(2001:2050,1:20) = 0;
bm(2051:2100,1:18) = 0;
bm(2101:2150,1:10) = 0;
bm(2151:2200,1:5) = 0;













bm(200:250,108:end) =0 ;
bm(251:300,104:end) = 0;
bm(301:350,100:end) = 0;
bm(351:400,98:end) = 0;
bm(401:451,94:end) = 0;
bm(451:500,89:end) = 0;
bm(501:550,87:end) = 0;
bm(551:600,85:end) = 0;
bm(601:650,84:end) = 0;
bm(651:700,81:end) = 0;
bm(700:751,79:end) = 0;
bm(751:800,77:end) =0 ;
bm(800:850,75:end) = 0;
bm(851:900,73:end) = 0;
bm(901:950,73:end) = 0;
bm(951:1000,73:end) = 0;
bm(1001:1050,73:end) = 0;
bm(1051:1100,75:end) = 0;
bm(1100:1150,77:end) = 0;
bm(1151:1200,79:end) = 0;


bm(1201:1250,81:end) = 0;
bm(1251:1300,83:end) = 0;
bm(1301:1350,85:end) = 0;
bm(1351:1400,86:end) = 0;
bm(1401:1450,88:end) = 0;
bm(1451:1500,90:end) = 0;
bm(1501:1550,92:end) = 0;
bm(1551:1600,94:end) = 0;
bm(1601:1650,96:end) = 0;
bm(1651:1700,98:end) = 0;
bm(1701:1750,100:end) = 0;
bm(1751:1800,102:end) = 0;
bm(1801:1900,104:end) = 0;
bm(1901:1950,106:end) = 0;
bm(1951:2000,108:end) = 0;
bm(2001:2050,110:end) = 0;
bm(2051:2100,115:end) = 0;
bm(2101:2150,120:end) = 0;
bm(2151:2200,125:end) = 0;

figure;
%subplot(132)
imagesc(bm);colormap gray
title('mask')


%% implementing mask to each transmit

for k = 1:128
        cd 'C:\Users\Valyria\Documents\Ultrasound\sumImage'
        %Files=dir('*.*');
        s = 'sumdata';
        num = num2str(k);
        FileNames = strcat(s,num);
%         Files=dir('*.*');
%         FileNames=Files(i+1).name;
        load(FileNames) 
        for j = 1:1:2353
    %for i = 1:1:128
            masking(j,:) = (sumsq(j,:).*bm(j,:));
            %disp(j)
            
    %end
        end
%         summask = sum(masking,2);
%         summasksq = squeeze(summask);
        %sumdat(:,k-2:k-2+127) = sumdat(:,k-2:k-2+127)+sumsq(:,:);
        s = 'sumMaskdata';
        num = num2str(k);
        filename = strcat(s,num);
        cd 'C:\Users\Valyria\Documents\Ultrasound\sumMask'
        save(filename,'masking')
        disp(k)
end



 %% adding mask data
 clc

sumMaskdat = zeros(2353,256);
% load('sumdata1.mat')
%summ = squeeze(sum(ImageData,2));
%sumdat(:,1:128) = sumsq(:,:);
% sumH = (20*log10(abs(hilbert(sumdat(2:2308,1:242))))); 
% sumH = -max(max(sumH))+sumH;
% 
% %subplot(121)
% imagesc(elementposmm,pixmm,sumH,[-50,0]);
% axis image
% colormap(gray)
% title('one beam')
for i = 1:128
        cd 'C:\Users\Valyria\Documents\Ultrasound\sumMask'
        s = 'sumMaskdata';
        num = num2str(i);
        FileNames = strcat(s,num);
%         Files=dir('*.*');
%         FileNames=Files(i+1).name;
        load(FileNames) 
        sumMaskdat(:,i:i+127) = sumMaskdat(:,i:i+127)+masking(:,:);
         
end
%    sumIm = sum(im_data,2);
%    sumIm2 = squeeze(sumIm);
ele_256posmm = (-255/2*es:es:255/2*es)*1e-3;
   im_scan_mask =(20*log10(abs(hilbert(sumMaskdat(2:2308,:))))); 
   im_scan_mask = -max(max(im_scan_mask))+im_scan_mask;
% 
     subplot(122)
    imagesc(ele_256posmm,pixmm,im_scan_mask,[-50,0]);
    axis image
    colormap(gray)
    title('masked data')

% for i = 1:1:2353
%     vals = elementnpos< round(elementnpos(i));
%     bwmask(i,vals) = 0;
% end



%subplot(121)

%elementposmm=(-127/2*es:es:127/2*es)*10^-3;
% masksum = sum(masking,2);
% mask_sq = squeeze(masksum);





%% wrong mask
first_half = mn(1:end,1:100)<=-40;
np = zeros(2307,256);
np = (first_half==1);
for i = 1:2307
    for j = 1:100
        if np(i,j) == 1
            bm(i,j) = 0;
        else
            bm(i,j) = 1;
        end
    end
end
second_half = mn(1:end,100:256)<=-40;
hj = (second_half==1);
gg = size(hj);
gk = gg(2);
for i = 1:2307
    for j = 1:gk
        if hj(i,j) == 1
            bm(i,100+j) = 0;
        else
            bm(i,100+j) = 1;
        end
    end
end

figure;
imagesc(bm)
colormap gray
for j = 1:1:2307
    for i = 1:1:256
            masking(j,i) = (sumdat(j+1,i).*bm(j,i));
            
            
    end
end

mask_sq2 = (20*log10(abs(hilbert(masking))));
bmode_mask = -max(max(mask_sq2))+mask_sq2;
figure;
imagesc(ele_256posmm,pixmm,bmode_mask,[-50,0])
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
            


