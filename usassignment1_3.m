clc;
clear all
close all

%%
cd('C:\Users\Valyria\Documents\Ultrasound')
load('pointTargetData.mat')
%%
% load('anecoicCystData.mat')
%%
no_lines = 128;
f0=veraStrct.frequencyMHz*1e6;%center frequency
fs = veraStrct.samplingRateMHz*1e6;
c = 1540;
data_pt = veraStrct.data(80:end,:,:);
depth = data_pt(:,1:128,65);%choosing lateral location at 65 - A-line 65
% index_channel = data_pt(1:end,2);
% lateral_location = data_pt(1:end,3);
pixelintomm = (0.5*c)/fs;
% min_depth = min(min(depth))*pixelintomm;
%depth_x = (0:size(depth,1)+min_depth);
pixnum = 1:1:length(depth);
pixmm = pixnum *pixelintomm;

% max_depth = max(max(depth))*pixelintomm;

imagesc(1:128,pixmm,depth,[-100 100])
colormap gray
xlabel('tranducer elements')
ylabel('depth(m)')



%point_image = [0 0 70]; %point to be imaged



x = 64.5;
dx = veraStrct.elementSpacingMM;%maybe center to center spacing too?

num_elements = 128;




test=1:length(data_pt(:,1,1));
es=veraStrct.elementSpacingMM;

elementposmm=(-127/2*es:es:127/2*es)*10^-3;

%do we create a new emit and receiver aperture
% emit_aperture = xdc_linear_array (n_elements, lambda/2, element_height, kerf, 1, 1,focus);
% xdc_impulse(emit_aperture,depth);
%% part 1
focal_point = [64.5,40*1e-3];
focal_p_1 = [64.5,40*1e-3];
data_num =  1:2353;
xi = linspace(1,2353,2392);

%data3 = interp1(1:1:2353,data_pt(:,:,:),xi,'linear');%linearly upsample the data
%upsample_data(:,:) = resample(data_pt(:,:),3,2);
depthmm = linspace(0.01,0.09,2353);
%calculating tau
for j = 1:1:num_elements
    data_fp = 64.5;
    td(j) = sqrt((elementposmm(j)^2+focal_p_1(2)^2))./c;
end

tdnew = td-max(td)*ones(size(td));

tdinpix = round(tdnew*fs*10);
data = data_pt(:,:,:);
output = zeros(size(data,1)*10-9,num_elements,num_elements);
%calculation of delayed data
for i = 1:num_elements
    line = squeeze(data(:,:,i));
    line_shift = zeros(size(data,1)*10-9,num_elements);
    tshift = line_shift;
    for j = 1:num_elements
        line_shift(:,j) = interp1(1:length(data),line(:,j),1:0.1:length(data));
        tshift(abs(tdinpix(j))+1:end,j) = line_shift(1:end+tdinpix(j),j);
    end
    output(:,:,i) = tshift;
end

imagesc(output(:,:,65),[-100,100]);colormap gray
% figure;
% 
% %have 3 dimensional data turn into 2 dimensional data
% 
% xlabel('element number')
% ylabel('pixels')
% title('part 1')

subplot (121)
sumnoInterp = sum(data_pt,2);
sumnoInterp = squeeze(sumnoInterp);
bmode_noInterp = (20*log10(abs(hilbert(sumnoInterp'))));
imagesc(elementposmm, depthmm,bmode_noInterp');
colormap(gray)
title('part 1 no delayed data')
axis image
% 
% subplot(212)
% sumInterp = sum(interpdata,2);
% sumInterp = squeeze(sumInterp);
% bmode = (20*log10(abs(hilbert(sumInterp'))));
% imagesc(bmode');
% colormap(gray)
% title('part 1 delayed data')

 subplot(122)
% %sumInterp = sum(interpdata,2);
sumInterp = sum(output,2);
sumInterp = squeeze(sumInterp);
bmode = (20*log10(abs(hilbert(sumInterp'))));
imagesc(elementposmm, depthmm,bmode');
colormap(gray)
title('part 1 delayed data')

%% part 2 focal point 1
focal_point_p1 = [64.5,10*1e-3];

for j = 1:1:num_elements
    td(j) = sqrt((elementposmm(j)^2+focal_point_p1(2)^2))./c;
end
%tdnew = td-min(td);

tdinpix = (td*fs);


for j = 1:1:128
    offset1 = min(td);
    td_f1(j) = sqrt((elementposmm(j)^2+focal_point_p1(2)^2))./c+(offset1/c);
end

td_f1_inpix=(td_f1-min(td_f1)).*fs;
mintdinpix = min(td)*fs;


%% part 2 focal point 2

focal_point2 = [64.5,30*1e-3];


test=1:length(data_pt(:,1,1));
es=veraStrct.elementSpacingMM;



%do we create a new emit and receiver aperture
% emit_aperture = xdc_linear_array (n_elements, lambda/2, element_height, kerf, 1, 1,focus);
% xdc_impulse(emit_aperture,depth);


for j = 1:1:num_elements
  
    td2(j) = sqrt((elementposmm(j)^2+focal_point2(2)^2))./c;
end

for j = 1:1:128
    offset2 = min(td2);
    td_f2(j) = sqrt((elementposmm(j)^2+focal_point2(2)^2))./c+(offset2/c);
end

td_f2_inpix=(td_f2-min(td_f2)).*fs;
mintd2inpix = min(td2)*fs;
% for i = 1:1:128
% for j = 521:1:1041
%  
%      td_f2(j-521+1,i) = sqrt((elementposmm(i)^2+focal_point2(2)^2))./c+(2*pixmm(j)/c);
% end
% 
% end
tdinpix2 = (td2*fs);



 

% for i = 1:num_elements
%     for j = 1:num_elements
%         interpdata2(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(tdinpix2(i))');
%     end
% end


%have 3 dimensional data turn into 2 dimensional data
% sumInterp2 = sum(interpdata2,2);
% sumInterp2 = squeeze(sumInterp2);
% bmode2 = (20*log10(abs(hilbert(sumInterp2'))));
% imagesc(bmode2');
% colormap(gray)
%% part 2 focal point 3
focal_point3 = [64.5,50*1e-3];


test=1:length(data_pt(:,1,1));
es=veraStrct.elementSpacingMM;



%do we create a new emit and receiver aperture
% emit_aperture = xdc_linear_array (n_elements, lambda/2, element_height, kerf, 1, 1,focus);
% xdc_impulse(emit_aperture,depth);


for j = 1:1:num_elements
  
    td3(j) = sqrt((elementposmm(j)^2+focal_point3(2)^2))./c;
end
for j = 1:1:128
    offset3 = min(td3);
    td_f3(j) = sqrt((elementposmm(j)^2+focal_point3(2)^2))./c+(offset3/c);
end


mintd3inpix = min(td3)*fs;
% for i = 1:1:128
% for j = 1042:1:1559
%      td_f3(j-1042+1,i) = sqrt((elementposmm(i)^2+focal_point3(2)^2))./c+(2*pixmm(j)/c);
% end
% end
tdinpix3 = (td3*fs);
td_f3_inpix=(td_f3-min(td_f3)).*fs;
% for i = 1:num_elements
%     for j = 1:num_elements
%         interpdata3(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(tdinpix3(i))');
%     end
% end


%have 3 dimensional data turn into 2 dimensional data
% sumInterp3 = sum(interpdata3,3);
% sumInterp3 = squeeze(sumInterp3);
% bmode3 = (20*log10(abs(hilbert(sumInterp3'))));
% imagesc(bmode3');
% colormap(gray)
%% part 2 focal point 4

focal_point4 = [64.5,70*1e-3];





%do we create a new emit and receiver aperture
% emit_aperture = xdc_linear_array (n_elements, lambda/2, element_height, kerf, 1, 1,focus);
% xdc_impulse(emit_aperture,depth);


for j = 1:1:num_elements
  
    td4(j) = sqrt((elementposmm(j)^2+focal_point4(2)^2))./c;
end
for j = 1:1:128
    offset4 = min(td4);
    td_f4(j) = sqrt((elementposmm(j)^2+focal_point4(2)^2))./c+(offset4/c);
end

mintd4inpix = min(td4)*fs;
% for i = 1:1:128
% for j = 1560:1:2079
%      td_f4(j-1560+1,i) = sqrt((elementposmm(i)^2+focal_point4(2)^2))./c+(2*pixmm(j)/c);
% end
% end
tdinpix4 = (td4*fs);


td_f4_inpix=(td_f4-min(td_f4)).*fs;

% for i = 1:num_elements
%     for j = 1:num_elements
%         interpdata4(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(tdinpix4(i))');
%     end
% end
% 
% 
% %have 3 dimensional data turn into 2 dimensional data
% sumInterp4 = sum(interpdata4,2);
% sumInterp4 = squeeze(sumInterp4);
% bmode4 = (20*log10(abs(hilbert(sumInterp4'))));
% imagesc(bmode4');
% colormap(gray)

%% part 2 focal point 5
focal_point5 = [64.5,85*1e-3];


test=1:length(data_pt(:,1,1));
es=veraStrct.elementSpacingMM;



%do we create a new emit and receiver aperture
% emit_aperture = xdc_linear_array (n_elements, lambda/2, element_height, kerf, 1, 1,focus);
% xdc_impulse(emit_aperture,depth);


for j = 1:1:num_elements
  
    td5(j) = sqrt((elementposmm(j)^2+focal_point5(2)^2))./c;
end
for j = 1:1:128
    offset5 = min(td5);
    td_f5(j) = sqrt((elementposmm(j)^2+focal_point5(2)^2))./c+(offset5/c);
end


mintd5inpix = min(td5)*fs;
% for i = 1:1:128
% for j = 2080:1:2353
%      td_f5(j-2080+1,i) = sqrt((elementposmm(i)^2+focal_point5(2)^2))./c+(2*pixmm(j)/c);
% end
% end
tdinpix5 = (td5*fs);
td_f5_inpix=(td_f5-min(td_f5)).*fs;
% for i = 1:num_elements
%     for j = 1:num_elements
%         interpdata5(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(tdinpix5(i))');
%     end
% end 
% 
% 
% %have 3 dimensional data turn into 2 dimensional data
% sumInterp5 = sum(interpdata5,2);
% sumInterp5 = squeeze(sumInterp5);
% bmode5 = (20*log10(abs(hilbert(sumInterp5'))));
% imagesc(bmode5');
% colormap(gray)

total_tdinpix = [td_f1_inpix;td_f2_inpix;td_f3_inpix;td_f4_inpix;td_f5_inpix];
% time_array = linspace(0*(fs/c),0.09*(fs/c),2353);

%this is only for region 1 thru 20


%%
original_timearray = ((1:1:2353)./fs);

%time array for each channel?! channel one foes from 1 until 2353 , channel
%2 goes from 1 until 2353
time_array = zeros(128,2353);

for i = 1:1:128
    for j = 1:1:2353
        time_array(i,j) = j;
    end
end
     

for i = 1:1:num_elements
    for j =1:1:num_elements
        interpdata1(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(td_f1_inpix(i)));
        interpdata2(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(td_f2_inpix(i)));
        interpdata3(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(td_f3_inpix(i)));
        interpdata4(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(td_f4_inpix(i)));
        interpdata5(:,i,j) = interp1((1:1:length(depth))',data_pt(:,i,j),(1:1:length(depth))'+round(td_f5_inpix(i)));
    end
end


 
%suminterp =[squeeze(sum(interpdata(1:round(mintdinpix),2)));squeeze(sum(interpdata2(round(mintd2inpix)+1:round(mintd3inpix),2)));squeeze(sum(interpdata3(round(mintd3inpix)+1:round(mintd4inpix),2)));squeeze(sum(interpdata4(round(mintd4inpix)+1:round(mintd5inpix),2)));squeeze(sum(interpdata5(round(mintd5inpix)+1:2353),2))];
suminterp1 = sum(interpdata1,2);
suminterp1 = squeeze(suminterp1);

suminterp2 = sum(interpdata2,2);
suminterp2 = squeeze(suminterp2);

suminterp3 = sum(interpdata3,2);
suminterp3 = squeeze(suminterp3);

suminterp4 = sum(interpdata4,2);
suminterp4 = squeeze(suminterp4);

suminterp5 = sum(interpdata5,2);
suminterp5 = squeeze(suminterp5);

%% part 2 calculation
part1 = suminterp1(1:1:520,:);
part2 = suminterp1(520:1:1041,:);
part3 = suminterp1(1041:1:1559,:);
part4 = suminterp1(1559:1:2079,:);
part5 = suminterp1(2079:1:2353,:);
sumInterpTotal=[part1;part2;part3;part4];



%% part 3
% original_timearray = ((1:1:2353)./fs);
% time_array = 1:(1/fs):original_timearray(2353);
focal_point_es = linspace(0,90*1e-3,2353);
for i = 1:1:length(focal_point_es)
    for j = 1:1:128
        fp = focal_point_es(i);
        td_es(i,j) = sqrt((elementposmm(j)^2+fp^2))./c+(fp/c);
        
    end
end

for i = 1:1:128
    td_es(i,:) = td_es(i,:);
end
td_es_pix = td_es*fs;
td_sub = (1:1:length(depth))'+td_es_pix;


for i = 1:1:length(depth)
    for k = 1:1:128
        delays = (1:1:length(depth))'+round(td_es_pix(i,k));
    end
end

time2 =time_array';
%delays2 = time_array+round(td_es_pix);

for i = 1:1:num_elements
    for j = 1:1:num_elements
         interpdata_es(:,i,j) = interp1(time2(:,i),data_pt(:,i,j),round(td_es_pix(:,i)));%try adding time2 to it
    end
end

suminterp_es = sum(interpdata_es,2);
suminterp_es2 = squeeze(suminterp_es);

%% plotting part 2 and part 3

subplot(121)
% sumInterp_part2 = sum(sumInterpTotal,2);
% sumInterp_part2 = squeeze(sumInterp_part2);
bmode_part2 = (20*log10(abs(hilbert(sumInterpTotal'))));
imagesc(elementposmm, depthmm,bmode_part2');
colormap(gray)
title('part 2')
xlabel('depth')
ylabel('element position in mm')
axis image


subplot(122)
% sumInterp_part3 = sum(suminterp_es2,2);
% sumInterp_part3 = squeeze(sumInterp_part3);
bmode_part3 = (20*log10(abs(hilbert(suminterp_es2'))));
imagesc(elementposmm, depthmm,bmode_part3');
colormap(gray)
title('part 3')
axis image

%% part 5 
focal_point_p5 = linspace(0,90*1e-3,2353);
data_pt2 = zeros(2353,128,64);

pitch = (veraStrct.elementSpacingMM)*1e-3;
%xf = (-8:1:8)*(pitch/2);
xf = [-0.5*pitch,0.5*pitch,-pitch,pitch,-1.5*pitch,1.5*pitch,-2*pitch,2*pitch,-2.5*pitch,2.5*pitch,-3*pitch,3*pitch,-3.5*pitch,3.5*pitch,-4*pitch,4*pitch];
xf1 = -pitch/2;
xf2 = pitch/2;
for j = 1:1:128
    for i = 1:1:length(focal_point_p5)
        
            delay1(i,j) = real((sqrt((elementposmm(j)^2-xf(1)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay2(i,j) = real((sqrt((elementposmm(j)^2-xf(2)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
            delay3(i,j) = real((sqrt((elementposmm(j)^2-xf(3)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay4(i,j) = real((sqrt((elementposmm(j)^2-xf(4)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
        
            delay5(i,j) = real((sqrt((elementposmm(j)^2-xf(5)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay6(i,j) = real((sqrt((elementposmm(j)^2-xf(6)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
            delay7(i,j) = real((sqrt((elementposmm(j)^2-xf(7)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay8(i,j) = real((sqrt((elementposmm(j)^2-xf(8)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
            delay9(i,j) = real((sqrt((elementposmm(j)^2-xf(9)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay10(i,j) = real((sqrt((elementposmm(j)^2-xf(10)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
            delay11(i,j) = real((sqrt((elementposmm(j)^2-xf(11)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay12(i,j) = real((sqrt((elementposmm(j)^2-xf(12)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
            delay13(i,j) = real((sqrt((elementposmm(j)^2-xf(13)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay14(i,j) = real((sqrt((elementposmm(j)^2-xf(14)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
            delay15(i,j) = real((sqrt((elementposmm(j)^2-xf(15)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            delay16(i,j) = real((sqrt((elementposmm(j)^2-xf(16)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
%             delay17(i,j) = real((sqrt((elementposmm(j)^2-xf(8)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay18(i,j) = real((sqrt((elementposmm(j)^2-xf(26)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             
%             delay19(i,j) = real((sqrt((elementposmm(j)^2-xf(7)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay20(i,j) = real((sqrt((elementposmm(j)^2-xf(27)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             
%             delay21(i,j) = real((sqrt((elementposmm(j)^2-xf(6)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay22(i,j) = real((sqrt((elementposmm(j)^2-xf(28)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             
%             delay23(i,j) = real((sqrt((elementposmm(j)^2-xf(5)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay24(i,j) = real((sqrt((elementposmm(j)^2-xf(29)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             
%             delay25(i,j) = real((sqrt((elementposmm(j)^2-xf(4)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay26(i,j) = real((sqrt((elementposmm(j)^2-xf(30)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             
%             delay27(i,j) = real((sqrt((elementposmm(j)^2-xf(3)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay28(i,j) = real((sqrt((elementposmm(j)^2-xf(31)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             
%             delay29(i,j) = real((sqrt((elementposmm(j)^2-xf(2)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay30(i,j) = real((sqrt((elementposmm(j)^2-xf(32)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             
%             delay31(i,j) = real((sqrt((elementposmm(j)^2-xf(1)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
%             delay32(i,j) = real((sqrt((elementposmm(j)^2-xf(33)^2)+focal_point_p5(i)^2)/c)+(focal_point_p5(i)/c))*fs;
            
           
    end
end
% delay1 = real(delay1)*fs;
% delay2 = real(delay2)*fs;
T = time_array';
%% for n = 2
pl = 2;
data = data_pt(:,:,pl/2:pl:end-pl/2);
x = size(data);
m = 1:1:2353;
for i = 1:x(3)
    DataForALine = squeeze(data(:,:,i));
        for j = 1:128
            SingleChannelData = DataForALine(:,j);
            ChDataDelay1 = interp1(T(:,j), SingleChannelData, delay1(:,j));
            ImageData(:,j,pl*i-1) = ChDataDelay1;

            ChDataDelay2 = interp1(T(:,j), SingleChannelData, delay2(:,j));
            ImageData(:,j,pl*i-0) = ChDataDelay2;
        end
end
% subplot(141)
figure;
    sumNewImage = sum(ImageData,2);
    sumNewImage = squeeze(sumNewImage);
    bmode_es5 = (20*log10(abs(hilbert(sumNewImage'))));
    bmode_es5 = -max(max(bmode_es5))+bmode_es5;
    imagesc(1:128,pixmm,bmode_es5',[-100 0]);
    colormap(gray)
    xlabel('number of elements')
    ylabel('depth')
    title('part 5: N =2 point target')
    
    



%% plots for all parts
% figure;
% subplot(131)
% imagesc(1:128,pixmm,interpdata(:,:,65),[-100 100])
% colormap gray
% xlabel('tranducer elements')
% ylabel('depth(m)')
% title('part 1')
% 
% subplot(132)
% bmode6 = (60*log10(abs(hilbert(sumInterpTotal'))));
% bmode6 = -max(max(bmode6))+bmode6;
% imagesc(1:128,pixmm,bmode6');
% colormap(gray)
% 
% xlabel('number of elements')
% ylabel('depth')
% title('part 2')
% %pbaspect([1 1 1])
% 
% subplot(133)
% bmode_es = (20*log10(abs(hilbert(suminterp_es2'))));
% bmode_es2 = -max(max(bmode_es))+bmode_es;
% imagesc(1:128,pixmm,bmode_es2',[-50 0]);
% colormap(gray)
% xlabel('number of elements')
% ylabel('depth')
% title('part 3')


% for i = 1:1:128
%     pitchposm(i) = elementposmm(i)-(pitch/2);
% end
% % distance_pitch = (elementposmm(1)-(pitch/2)):pitch:(elementposmm(128)-(pitch/2));
% % distance_pitch = distance_pitch.*1e-3;
% 
% 
%     %for k =2:1:length(distance_pitch)
% 
% tau1inpix = real(tau1)*fs;
% tau2inpix = real(tau2)*fs;
% % parallel beams 4
% 
% for i = 1:1:128
%     for j = 1:1:4
%         interpt1data(:,i,j) = interp1(time2(:,i)', data_pt2(:,i,j),tau1inpix(:,i));
%     end
% end
% 
% interp1datasum = sum(interpt1data,2);
% interp1datasum = squeeze(interpt1datasum);
% 
% bmode_par = (60*log10(abs(hilbert(interp1datasum'))));
% imagesc(1:128,pixmm,bmode_par');
% colormap(gray)

%pbaspect([1 1 1])



% figure;
% 
% subplot(231)
% bmode7 = (20*log10(abs(hilbert(suminterp1'))));
% 
% 
% imagesc(1:128,pixmm,bmode7');
% colormap(gray)
% xlabel('number of elements')
% ylabel('depth')
% title('image at first focal depth')
% 
% subplot(232)
% bmode8 = (20*log10(abs(hilbert(suminterp2'))));
% imagesc(1:128,pixmm,bmode8');
% colormap(gray)
% % yaxis = ((1:1:2353)*(c/fs));
% % axis([1 128 yaxis(1) yaxis(2353)])
% xlabel('number of elements')
% ylabel('depth')
% title('image at second focal depth')
% 
% 
% subplot(233)
% bmode9 = (20*log10(abs(hilbert(suminterp3'))));
% imagesc(1:128,pixmm,bmode9');
% colormap(gray)
% % yaxis = ((1:1:2353)*(c/fs));
% % axis([1 128 yaxis(1) yaxis(2353)])
% xlabel('number of elements')
% ylabel('depth')
% title('image at third focal depth')
% 
% subplot(234)
% bmode10 = (20*log10(abs(hilbert(suminterp4'))));
% imagesc(1:128,pixmm,bmode10');
% colormap(gray)
% % yaxis = ((1:1:2353)*(c/fs));
% % axis([1 128 yaxis(1) yaxis(2353)])
% xlabel('number of elements')
% ylabel('depth')
% title('image at fourth focal depth')
% 
% subplot(235)
% bmode11 = (20*log10(abs(hilbert(suminterp5'))));
% imagesc(1:128,pixmm,bmode11');
% colormap(gray)
% % yaxis = ((1:1:2353)*(c/fs));
% % axis([1 128 yaxis(1) yaxis(2353)])
% xlabel('number of elements')
% ylabel('depth')
% title('image at fifth focal depth')
% 
% subplot(236)
% bmode6 = (20*log10(abs(hilbert(sumInterpTotal'))));
% imagesc(1:128,pixmm,bmode6');
% colormap(gray)
% % yaxis = ((1:1:2353)*(c/fs));
% % axis([1 128 yaxis(1) yaxis(2353)])
% xlabel('number of elements')
% ylabel('depth')
% title('part 2')



% 
% bmode5 = (20*log10(abs(hilbert(suminterp'))));
% imagesc(bmode5');
% colormap(gray)
%  
% for j = 1:1:128
%     for i = 1:1:2353
%         time_array(j,i) = i;
%     end
% end
%  for i = 1:num_elements
%   for j= 1:1:2353  
%      interpfordata(j) = (1:1:length(time_array))'+round(total_tdinpix(i,j));
%      %interpdata_total(:,i) = interp1((1:1:length(time_array))',data_pt(:,i),interpfordata);
%   end
%  end
%  %interpolated data for 1 focal point
%  point = total_tdinpix(1,:);
%  for i = 1:num_elements
%     for j = 1:num_elements
%         interp_point(:,i,j) = interp1((1:1:length(time_array))',data_pt(:,i,j),(1:1:length(time_array))'+total_tdinpix(1,i));
%     end
%  end
% 
% suminterp = sum(interp_point,2);
% suminterp = squeeze(suminterp);
% 
% bmode5 = (20*log10(abs(hilbert(inter_channel128'))));
% imagesc(bmode5');
% colormap(gray)
%  
%  
% inter_channel1 = interp1((1:1:length(time_array))',data_pt(:,1),(1:1:length(time_array))'+round(total_tdinpix(1,:)));
% inter_channel2 = interp1((1:1:length(time_array))',data_pt(:,2),(1:1:length(time_array))'+round(total_tdinpix(2,:)));
% inter_channel128 = interp1((1:1:length(time_array))',data_pt(:,128),(1:1:length(time_array))'+round(total_tdinpix(128,:)));
% 
% sumInterptotal = sum(interpdata_total,2);
% sumInterptotal = squeeze(sumInterptotal);
% 
% bmode5 = (20*log10(abs(hilbert(inter_channel128'))));
% imagesc(bmode5');
% colormap(gray)
