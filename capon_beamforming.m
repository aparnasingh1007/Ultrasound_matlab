clear all; close all;
fs = 10e6; %sampling frequency
ts = 1/fs; 
T = 0.0001; %signal duration
t = 0:ts:T-ts;
Ns = length(t); %number of samples
fc = 1e6; %carrier freqiency
lambda = 3*10^8/fc; %wavelength
M = 10;%number of antenna elements in Unliform Linear Array
d = 0.5*lambda;%distance between neigbouring antennas
%------------------------------------------------------------------
s1 = exp(1i*2*pi*fc*t); %signal arriving to the array phase center (first antenna) - recieve signal
load('channel_data.mat')
%%
% for i = 1:5
%     cd 'C:\Users\Valyria\Documents\Ultrasound\interpdata'
%     s = 'interpdata';
%     num = num2str(i);
%    filename = strcat(s,num);
%    load(filename)
%    recieve_data(:,:,i) = ImageData(:,:);%receive data after synthetic apeture
% end
cd 'C:\Users\Valyria\Documents\Ultrasound\interpdata'
load('interpdata1.mat')
recieve_data1 = ImageData(:,:,1);
%%
%SINGLE ARRIVING SIGNAL - WORKS FINE BOTH FOR CONVENTIONAL AND CAPON
%BEAMFORMERS
%teta = [40]/180*pi; %direction of arrival in degrees 
%TWO ARRIVING SIGNALS - CONVENTIONAL BF - OK, CAPON BF - INVALID SPECTRUM!!!
%%
teta = [0, 40]/180*pi; %directions of arrivals in degrees 
amp = [1 1];% amplitudes of signals
delta_fi = -2*pi*d/lambda*sin(teta); %relative phase delay between signals received through neogbouring elements
for m=1:M;
    aU(:,m) = exp(1i*((m-1)*delta_fi)); %steering vector
end
%%
%x = zeros(M,Ns);
% for k=1:length(teta);
%      x = x + amp(k)*aU(k,:).'*exp(1i*randn(size(t))); %data at the outputs of antenna elements - trasnmit data 
% end
%image data
x = chanOutStore(:,:,1);
SNR = 10; %signal-to-noise ratio in decibels
sz = (2^(-0.5))*sqrt(10^(-SNR/10))*(randn(M,Ns)+1i*randn(M,Ns)); %complex Gaussian noise
%%
Px = var(aU(k,1).'*s1*exp(1i*randn())); %power of signal
x = x + sz; %adding noise
Psz = var(sz(1,:)); %noise variance
Rxx = x*x'/Ns; %data covariance matrix
SNR_true = 10*log10(Px/Psz) %true SNR in dB
iRxx = inv(Rxx); %inverse of covariance matrix
teta_v = -pi/2:pi/180:pi/2; %range of scanned directions in radians
teta_v_deg = teta_v/pi*180; %same as above but in degrees
Nteta = length(teta_v); %number of scanned directions
%%
for k=1:Nteta; %scanning loop
    %k
    teta = teta_v(k); %current scannig direction
    delta_fi = -2*pi*d/lambda*sin(teta); %relateive phase delay 
    for m=1:M;
        aT(m) = exp(1i*((m-1)*delta_fi)); %steering vector
    end
      a=aT.'; %transpose of steering vector
      wbf = a; %weight vector for conventional beamforming
      Pbf(k) = wbf'*Rxx*wbf; %spatial power spectrum of conventional beamformer 
      wMVDR = (iRxx*a)/(a'*iRxx*a); %weight vector for Caponbeamforming
      PMVDR(k) = wMVDR'*Rxx*wMVDR/1; %spatial power spectrum of Capon beamformer
  end
Pbf_dB = 10*log10(Pbf); %spectrum in log scale
PMVDR_dB = 10*log10(PMVDR);
figure
plot(teta_v_deg, Pbf_dB,'k','LineWidth',2);
axis([-90 90 0 20*log10(M)+5])
xlabel('\theta [\circ]');
ylabel('P(\theta)','rotation',0);
figure
plot(teta_v_deg, PMVDR_dB,'k','LineWidth',2);
xlabel('\theta [\circ]');
ylabel('P(\theta)','rotation',0);