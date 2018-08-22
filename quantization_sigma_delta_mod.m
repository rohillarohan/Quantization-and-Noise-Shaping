%% Quantization and Noise Shaping using Sigma-Delta Modulation Tehcnique
clc
clearvars
close all
%% Part I: Quantization of a Signal
fsI = 29927; % Sampling Frequency
tsI = 1/fsI; % Sampling Interval
tI1 = 0:tsI:5000*tsI-tsI;
xI1 = 0.25*sin(3221*tI1) +0.25*cos(8169*tI1); % Signal
DI1 = 7; % Round-Off Factor
disp('------------------------------- Part:I --------------------------------------')
%% Part I: 1
% 
xQntI1 = round(xI1*2^DI1)/2^DI1; % Quantized Signal

qntErI1 = xI1-xQntI1; % Quantization Error

VI1 = var(qntErI1); % Average Power

%% Part I: 2
% 
acrQntErI1 = xcorr(qntErI1)/length(tI1); % Auto Correlation of the Error Signal

figure;
plot((-(length(acrQntErI1)-1)/2:(length(acrQntErI1)-1)/2),acrQntErI1,'linewidth', 1.5)
title('Autocorrelation of the Initial Error Signal')
xlabel('Samples')
ylabel('Sample Size')
grid

disp('Average power calculated (D = 7 & N = 5000):');disp(VI1)
disp('Average power observed from the plot for the initial signal:');disp(acrQntErI1(5000))

% Theoretical Power:
theorQNPwrI1 = 1/(2^(2*DI1)*12);
disp('Theoretical Power:');disp(theorQNPwrI1)

% Why is theoretical Quantization Noise power different from the one we
% calculated?
% 
% Explanation: The reason why theoretical power is less than the one we
% calculated because MATLAB also quantizes the values although it is called
% infinite precision (64 bit in this case) it is actually not.

%% Part I: 3
figure;

% PDF histogram of the Error Signal Samples
hist(qntErI1)

title('PDF histogram of the Initial Error Signal Samples')

% Why the Distribution of the values is not exactly uniform?
% 
% Explanation: Here we should note that the the theoretical model of
% quantization noise assumes that its amplitude values follows the uniform
% probability density function. But actually the quantization error is not
% uniformlly distributed. This is because as we select the number of bits
% (D) for quantization, the fulscale voltage gets divided in 2^D
% quantization levels. Now the error between the ideal value and the
% quantized values ranges from -DELTA/2 to +DELTA/2 (DELTA = V_fullscale/128). The error could be
% anywhere between these limits. Therefore the distribution is not exactly
% uniform. But it is obviously near to uniform. It is pretty intuitive to
% think that if the defference between -DELTA/2 and +DELTA/2 decreases the 
% distribution gets more uniform in nature (because the difference between 
% two error values also decreases)

%% Part I: 4
tI2 = 0:tsI:30000*tsI-tsI;
xI2 = 0.25*sin(3221*tI2) +0.25*cos(8169*tI2); % Signal
DI2 = 7; % Round-Off Factor
% 
xQntI2 = round(xI2*2^DI2)/2^DI2; % Quantized Signal

qntErI2 = xI2-xQntI2; % Quantization Error

VI2 = var(qntErI2); % Average Power
% 
acrQntErI2 = xcorr(qntErI2)/length(tI2); % Auto Correlation of the Error Signal

figure;
plot((-(length(acrQntErI2)-1)/2:(length(acrQntErI2)-1)/2),acrQntErI2,'linewidth', 1.5)
title('Autocorrelation of the Error Signal When N = 30000')
xlabel('Samples')
ylabel('Sample Size')
grid

disp('Average power calculated when N = 30000, for the second signal:');disp(VI2)
disp('Average power observed from the plot for the second signal:');disp(acrQntErI2(30000))

figure;
hist(qntErI2)
title('PDF histogram When N = 30000, D = 7')

% Explain the changes.
% 
% Explanation: Here as we increased the number of samples the changes we
% noted are that the distribution of quantized noise got more uniformly
% distributed and the power of the Quantization noise was not effected much (5.2e-6 to 5.05e-6).
% The reason why this happened lies in the fact that we increased the
% number of samples and we did not increase the number of bits. One may ask
% what is the difference? The difference is that when we increase the number of
% samples the ratio of the difference between two samples of the
% quantization noise to the number of samples decreases. This makes the
% distribution looks more uniform. Note that the ratio decreases not the
% difference between two quantization noise samples. This difference is inversely 
% propotional to the number of bits used to quantiz the signal. As we did not increase
% the number of bits the difference did not reduce. This is in fact supported
% by our second observation i.e. the powerof the quantization noise, as we did not
% increase the number of bits used for quantization the power barely reduced. 

%% Part I: 5
tI3 = 0:tsI:5000*tsI-tsI;
xI3 = 0.25*sin(3221*tI3) +0.25*cos(8169*tI3); % Signal
DI3 = 12; % Round-Off Factor
% 
xQntI3 = round(xI3*2^DI3)/2^DI3; % Quantized Signal

qntErI3 = xI3-xQntI3; % Quantization Error

VI3 = var(qntErI3); % Average Power
% 
acrQntErI3 = xcorr(qntErI3)/length(tI3); % Auto Correlation of the Error Signal

figure;
plot((-(length(acrQntErI3)-1)/2:(length(acrQntErI3)-1)/2),acrQntErI3,'linewidth', 1.5)
title('Autocorrelation of the Error Signal When D = 12')
xlabel('Samples')
ylabel('Sample Size')
grid

disp('Average power calculated when D = 12 for the third signal:');disp(VI3)
disp('Average power observed from the plot for the third signal:');disp(acrQntErI3(5000))

figure;
hist(qntErI3) % PDF histogram of the Error Signal Samples
title('PDF histogram When N = 5000,  D = 12')

% Explain the changes in the values.
% 
% Explanation: Now from above two observations and their comprehensions we
% can say that the changes we see which are the distribution is more
% uniform and the poer of the quantization noise is reduced to a much lower value
% (5.2e-6 to 4.94e-9) were pretty much expected. As we increased the number
% of bits used for quantization (7 to 12) we increased the number of quantization
% leves from (128 to 4096). The same ful scale voltage is now divided into
% 4096 levels. This reduced the limt of the quantization error -DELTA/2 to +DELTA/2 (DELTA = V_fullscale/4096). 
% as the limit reduced the difference between the values two noise sample
% can take is also decreased.

%% Part II: Finding quantization error at the output of a first order digital filter

fsII1 = 29927;
tsII1 = 1/fsII1;
tII1 = 0:tsII1:5000*tsII1-tsII1;
xII1 = 0.1*sin(2*pi*785*tII1)+0.1*sin(2*pi*1031*tII1); % Ideal Input

disp('------------------------------- Part:II D = 8,a = 0.6--------------------------------------')

%% Part II: 1
aII1 = .6;
HII1 = tf([1 0],[1 -aII1],1);
[ninf1,fpeak1] = norm(HII1,inf);
kII1 = 1/ninf1;
disp('K when a = 0.6');disp(kII1)

%% Part II: 2
figure;
subplot(221)
plot((0:length(xII1)-1),xII1,'linewidth',2);
title('Ideal Input Signal (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xII1)-.1 max(xII1)+.1])
grid

%% Part II: 3
DII1 = 8;
xQntII1 = round(xII1*2^DII1)/2^DII1; % Quantized Input

% Plotting the Quantized Input
subplot(222)
plot((0:length(xQntII1)-1),xQntII1,'linewidth',2);
title('Quantized Input Signal (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xQntII1)-.1 max(xQntII1)+.1])
grid

% Finding the Quantized Output
yTempII1 = zeros(1,5000);
yQntII1 = zeros(1,5000); % Quantized Output
yTempII1(1) = xQntII1(1);
yQntII1(1) = xQntII1(1);
for n = 1:4999  
    yTempII1(n+1) = (aII1*yTempII1(n))+(kII1*xQntII1(n+1));
    yQntII1(n+1) = round(yTempII1(n+1)*2^DII1)/2^DII1;
end

% Plotting the Quantized Output
subplot(223)
plot((0:length(yQntII1)-1),yQntII1,'linewidth',2);
title('Quantized Output Signal (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yQntII1)-.1 max(yQntII1)+.1])
grid

%% Part II: 4
% Calculating the Ideal Output
yIdlII1 = zeros(1,5000); 
yIdlII1(1) = kII1*xII1(1);
for n = 1:4999
    yIdlII1(n+1) = aII1*yIdlII1(n)+(kII1*xII1(n+1));
end

% Plotting the Ideal Output
subplot(224)
plot((0:length(yIdlII1)-1),yIdlII1,'linewidth',2);
title('Ideal Output Signal (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yIdlII1)-.1 max(yIdlII1)+.1])
grid

%% Part II: 5
yQntErII1 = yIdlII1-yQntII1; % Error between the outputs (Quantization Noise)
VErII1 = var(yQntErII1); % Power of the Quantization Noise
VQntSigII1 = var(yQntII1); % Power of the Quantized Output Signal
VIdlSigII1 = var(yIdlII1); % Power of the Ideal Output SignalSignal

% SnRII1 = 20*log(VQntSigII1/VErII1);
SnRII1a = 10*log10(VIdlSigII1/VErII1);
disp('Ideal Signal to Noise ratio (D = 8, a = 0.6):');disp(SnRII1a)

SnRII1b = 10*log10(VQntSigII1/VErII1);
disp('Quantization Signal to Noise ratio (D = 8, a = 0.6):');disp(SnRII1b)

%% Part II: 6
% Theoretical SNR Calculation
theorQNPwrII1 = 1/(2^(2*DII1)*12);
theorSnrII1 = 10*log10(VIdlSigII1/theorQNPwrII1);
disp('Theoretical: Ideal Signal to Noise ratio (D = 8, a = 0.6):');disp(theorSnrII1)
%% Part II: 7
% Plotting the Frequency Response of the Quantized Output

figure;
frqRspQntOPII1 = 20*log10(fftshift(abs(fft(yQntII1,4096))));
plot((-29927/2:29927/4096:29927/2 - (29927/4096)),frqRspQntOPII1);
title('Frequency Response of the Quantized Output (D = 8, a = 0.6)')
xlabel('Frequency (Hz)')
ylabel('Log Magnitude')
grid

% Explain why the level of the noise is higher at lower frequencies?
%
% Explanation: In case of the system with no noise shaping the noise
% transfer function is given by the equation [1/(1-a*(z^-1))] when we take
% the frequency transform (z = e^jw) we find that at lower frequencies the
% noise transfer function has higher values than at higher frequencies.
% This gives the higher noise at lower frequencies.

%% Part II: 8
disp('-------------------------------Part II: D = 12, a = 0.6--------------------------------------')
fsII2 = 29927;
tsII2 = 1/fsII2;
tII2 = 0:tsII2:5000*tsII2-tsII2;
xII2 = (0.1*sin(2*pi*785*tII2))+(0.1*sin(2*pi*1031*tII2)); % Ideal Input

aII2 = .6;
HII2 = tf([1 0],[1 -aII2],1);
[ninf2,fpeak2] = norm(HII2,inf);
kII2 = 1/ninf2;
disp('K when a = 0.6');disp(kII2)

% Plotting the Ideal Input
figure;
subplot(221)
plot((0:length(xII2)-1),xII2,'linewidth',2);
title('Ideal Input Signal (D = 12, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xII2)-.1 max(xII2)+.1])
grid

DII2 = 12;
xQntII2 = round(xII2*2^DII2)/2^DII2; % Quantized Input

% Plotting the Quantized Input
subplot(222)
plot((0:length(xQntII2)-1),xQntII2,'linewidth',2);
title('Quantized Input Signal (D = 12, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xQntII2)-.1 max(xQntII2)+.1])
grid

% Calculating Quantized Output
yTempII2 = zeros(1,5000);
yQntII2 = zeros(1,5000); % Quantized Output
yTempII2(1) = xQntII2(1);
for n = 1:4999  
    yTempII2(n+1) = (aII2*yTempII2(n))+(kII2*xQntII2(n+1));
    yQntII2 = round(yTempII2*2^DII2)/2^DII2;
end

% Plotting Quantized Output
subplot(223)
plot((0:length(yQntII2)-1),yQntII2,'linewidth',2);
title('Quantized Output Signal (D = 12, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yQntII2)-.1 max(yQntII2)+.1])
grid

% Calculating Ideal Output
yIdlII2 = zeros(1,5000); 
yIdlII2(1) = kII2*xII2(1);
for n = 1:4999
    yIdlII2(n+1) = (aII2*yIdlII2(n))+(kII2*xII2(n+1));
end

% Plotting the Ideal Output
subplot(224)
plot((0:length(yIdlII2)-1),yIdlII2,'linewidth',2);
title('Ideal Output Signal (D = 12, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yIdlII2)-.1 max(yIdlII2)+.1])
grid

yQntErII2 = yIdlII2-yQntII2; % Error between the outputs (Quantization Noise)
VErII2 = var(yQntErII2); % Power of the Quantization Noise
VQntSigII2 = var(yQntII2); % Power of the Quantized Output Signal
VIdlSigII2 = var(yIdlII2); % Power of the Ideal Output SignalSignal

SnRII2a = 10*log10(VIdlSigII2/VErII2);
disp('Ideal Signal to Noise ratio (D = 12, a = 0.6):');disp(SnRII2a)

SnRII2b = 10*log10(VQntSigII2/VErII2);
disp('Quantization Signal to Noise ratio (D = 12, a = 0.6):');disp(SnRII2b)

% Theoretical SNR Calculation
theorQNPwrII2 = 1/(2^(2*DII2)*12);
theorSnrII2 = 10*log10(VIdlSigII2/theorQNPwrII2);
disp('Theoretical: Ideal Signal to Noise ratio (D = 12, a = 0.6):');disp(theorSnrII2)

% Explain the differences between the SQNRs and the spectrums.
% 
% Explanation: When we increased the number of bits used to quantize the
% signal from 8 to 12 we increased the number of quantization levels in the
% fullscale voltage from 256 to 4096. This gave our quantization more
% precision as explained in first part. The limit of the quantization error
% got reduced to -DELTA/2 to +DELTA/2 (DELTA = V_fullscale/4096). This
% reduced the power of the Quantization noise and increased the SQNR
% significantly. This can be seen in the spectrum, Low noise is shown by
% lower spikes in the spectrum all over the frequency range.

% Plotting the Frequency Response Quantized Output
figure;
frqRspQntOPII2 = 20*log10(fftshift(abs(fft(yQntII2,4096))));
plot((-29927/2:29927/4096:29927/2 - (29927/4096)),frqRspQntOPII2);
title('Frequency Response of the Quantized Output (D = 12, a = 0.6)')
xlabel('Frequency (Hz)')
ylabel('Log Magnitude')
axis([-1.5e4 1.5e4 -40 60])
grid

%% Part II: 9
disp('-------------------------------Part II: D = 8,a = .95--------------------------------------')
fsII3 = 29927;
tsII3 = 1/fsII3;
tII3 = 0:tsII3:5000*tsII3-tsII3;
xII3 = 0.1*sin(2*pi*785*tII3)+0.1*sin(2*pi*1031*tII3); % Ideal Input

aII3 = .95;
HII3 = tf([1 0],[1 -aII3],1);
[ninf3,fpeak3] = norm(HII3,inf);
kII3 = 1/ninf3;
disp('K when a = 0.95');disp(kII3)

% Plotting Ideal Input
figure;
subplot(221)
plot((0:length(xII3)-1),xII3,'linewidth',2);
title('Ideal Input Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xII3)-.1 max(xII3)+.1])
grid

DII3 = 8;
xQntII3 = round(xII3*2^DII3)/2^DII3; % Quantized Input

% Plotting Quantized Input
subplot(222)
plot((0:length(xQntII3)-1),xQntII3,'linewidth',2);
title('Quantized Input Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xQntII3)-.1 max(xQntII3)+.1])
grid

% Calculating Quantized Output
yTempII3 = zeros(1,5000);
yQntII3 = zeros(1,5000); 
yTempII3(1) = kII3*xQntII3(1);
yQntII3(1) = kII3*xQntII3(1);
for n = 1:4999  
    yTempII3(n+1) = (aII3*yTempII3(n))+(kII3*xQntII3(n+1));
    yQntII3(n+1) = round(yTempII3(n+1)*2^DII3)/2^DII3;
end

% Plotting Quantized Output
subplot(223)
plot((0:length(yQntII3)-1),yQntII3,'linewidth',2);
title('Quantized Output Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yQntII3)-.1 max(yQntII3)+.1])
grid

% Calculating Ideal Output
yIdlII3 = zeros(1,5000); 
yIdlII3(1) = kII3*xII3(1);
for n = 1:4999
    yIdlII3(n+1) = (aII3*yIdlII3(n))+(kII3*xII3(n+1));
end

% Plotting Ideal Output
subplot(224)
plot((0:length(yIdlII3)-1),yIdlII3,'linewidth',2);
title('Ideal Output Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yIdlII3)-.1 max(yIdlII3)+.1])
grid

yQntErII3 = yIdlII3-yQntII3; % Error between the outputs (Quantization Noise)
VErII3 = var(yQntErII3); % Power of the Quantization Noise
VQntSigII3 = var(yQntII3); % Power of the Quantized Output Signal
VIdlSigII3 = var(yIdlII3); % Power of the Ideal Output SignalSignal

SnRII3a = 10*log10(VIdlSigII3/VErII3);
disp('Ideal Signal to Noise ratio (D = 8, a = 0.95):');disp(SnRII3a)

SnRII3b = 10*log10(VQntSigII3/VErII3);
disp('Quantization Signal to Noise ratio (D = 8, a = 0.95):');disp(SnRII3b)

% Theoretical SNR Calculation
theorQNPwrII3 = 1/(2^(2*DII3)*12);
theorSnrII3 = 10*log10(VIdlSigII3/theorQNPwrII3);
disp('Theoretical: Ideal Signal to Noise ratio (D = 8, a = 0.95):');disp(theorSnrII3)

% Why do we get worse SNR values compare to the case where the pole was 0.6?
% 
% Explaination: The reason why we got the worse SQNR as compared to the
% case when the pole was at 0.6 lies in the fact that for the given first
% order system the variance gain is inversely proportional the square of the value of 
% the pole [ 1/(1-(a^2)) ]. As close the pole gets to the unit circle (in the given case 0.95)
% we get a very high variance gain and thus a very high noise at the output.  

% Plotting the Frequency Response of Quantized Output
figure;
frqRspQntOPII3 = 20*log10(fftshift(abs(fft(yQntII3,4096))));
plot((-29927/2:29927/4096:29927/2 - (29927/4096)),frqRspQntOPII3);
title('Frequency Response of the Quantized Output (D = 8, a = 0.95)')
xlabel('Frequency (Hz)')
ylabel('Log Magnitude')
axis([-1.5e4 1.5e4 -50 50])
grid

%% Part III: Using noise shaping technique to increase the SNR
fsIII1 = 29927;
tsIII1 = 1/fsIII1;
tIII1 = 0:tsIII1:5000*tsIII1-tsIII1;
xIII1 = 0.1*sin(2*pi*785*tIII1)+0.1*sin(2*pi*1031*tIII1); % Ideal Input

disp('-----------Noise Shaping FilterPart:III D = 8,a = 0.6--------------')

%% Part III: 1:1
aIII1 = .6;
HIII1 = tf([2 -1],[1 -aIII1],1);
[ninf2,fpeak2] = norm(HIII1,inf);
kIII1 = 1/ninf2;
disp('K when a = 0.6');disp(kIII1)

%% Part III: 1:2
% Plotting Ideal Input
figure;
subplot(221)
plot((0:length(xIII1)-1),xIII1,'linewidth',2);
title('Ideal Input Signal (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xIII1)-.1 max(xIII1)+.1])
grid

%% Part III: 1:3
DIII1 = 8;
xQntIII1 = round(xIII1*2^DIII1)/2^DIII1; % Quantized Input

% Plotting the Quantized Input
subplot(222)
plot((0:length(xQntIII1)-1),xQntIII1,'linewidth',2);
title('Quantized Input Signal (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xQntIII1)-.1 max(xQntIII1)+.1])
grid

% Calculating Quantized Output
yTemp1III1 = zeros(1,5000);
yQntIII1 = zeros(1,5000); 
yTemp1III1(1) = round((kIII1*xQntIII1(1))*2^(2*DIII1))/2^(2*DIII1);
yQntIII1(1) = round((kIII1*xQntIII1(1))*2^(2*DIII1))/2^(2*DIII1);
yQntIII1(1) = round(yQntIII1(1)*2^(DIII1))/2^(DIII1);
err = zeros(1,5000);

for n = 1:4999 
    err(n) = round((-1*err(n)*2^(2*DIII1)))/(2^(2*DIII1));
    yTemp1III1(n+1) = (kIII1*xQntIII1(n+1))+(aIII1*yQntIII1(n))+err(n);
    yTemp1III1(n+1) = round(yTemp1III1(n+1)*2^(2*DIII1))/(2^(2*DIII1));
    yQntIII1(n+1) = round(yTemp1III1(n+1)*2^DIII1)/(2^DIII1);
    err(n+1) = round((yQntIII1(n+1)-yTemp1III1(n+1))*2^(2*DIII1))/(2^(2*DIII1));
end

% Plotting the Quantized Output
subplot(223)
plot((0:length(yQntIII1)-1),yQntIII1,'linewidth',2);
title('Quantized Output Signal for Noise Shaping Filter (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yQntIII1)-.1 max(yQntIII1)+.1])
grid

%% Part III: 1:4
% Calculating Ideal Output 
% As the system transfer function of the noise shaping filter is the same
% as the regular feedback system,
% We can take the ideal system response of the previous filter.
yIdlIII1 = zeros(1,5000);
yIdlIII1(1) = kIII1*xIII1(1);

for n = 1:4999
    yIdlIII1(n+1) = kIII1*xIII1(n+1)+(aIII1*yIdlIII1(n));
end

% Plotting the Ideal Output
subplot(224)
plot((0:length(yIdlIII1)-1),yIdlIII1,'linewidth',2);
title('Ideal Output Signal for Noise Shaping Filter (D = 8, a = 0.6)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yIdlIII1)-.1 max(yIdlIII1)+.1])
grid

%% Part III: 1:5
yQntErIII1 = yIdlIII1-yQntIII1; % Error between the outputs (Quantization Noise)
VErIII1 = var(yQntErIII1); % Power of the Quantization Noise
VQntSigIII1 = var(yQntIII1); % Power of the Quantized Output Signal
VIdlSigIII1 = var(yIdlIII1); % Power of the Ideal Output SignalSignal

SnRIII1a = 10*log10(VIdlSigIII1/VErIII1);
disp('Ideal Signal to Noise ratio (D = 8, a = 0.6):');disp(SnRIII1a)

SnRIII1b = 10*log10(VQntSigIII1/VErIII1);
disp('Quantization Signal to Noise ratio (D = 8, a = 0.6):');disp(SnRIII1b)

%% Part III: 1:6
% Theoretical SNR Calculation
theorQNPwrIII1 = 1/(2^(2*DIII1)*12);
theorSnrIII1 = 10*log10(VIdlSigIII1/theorQNPwrIII1);
disp('Theoretical: Ideal Signal to Noise ratio (Noise Shaping Filter) (D = 8, a = 0.6):');disp(theorSnrIII1)

%% Part III: 1:7
% Plotting the Frequency Responce of the Noise Shaping Filter
figure;
frqRspQntOPIII1 = 20*log10(fftshift(abs(fft(yQntIII1,4096))));
plot((-29927/2:29927/4096:29927/2 - (29927/4096)),frqRspQntOPIII1);
title('Frequency Response of the Quantized Output (D = 8, a = 0.6)')
xlabel('Frequency (Hz)')
ylabel('Log Magnitude')
grid

%% Part III: 3
disp('-----------Noise Shaping FilterPart:III D = 8,a = 0.95--------------')
aIII2 = .95;
HIII2 = tf([2 -1],[1 -aIII2],1);
[ninf3,fpeak3] = norm(HIII2,inf);
kIII2 = 1/ninf3;
disp('K when a = 0.95');disp(kIII2)

xIII2 = xIII1;

DIII2 = 8;
xQntIII2 = round(xIII2*2^DIII2)/2^DIII2; % Quantized Input

% Plotting the Ideal Input
figure;
subplot(221)
plot((0:length(xIII2)-1),xIII2,'linewidth',2);
title('Ideal Input Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xIII2)-.1 max(xIII2)+.1])
grid

% Plotting the Quantized Input
subplot(222)
plot((0:length(xQntIII2)-1),xQntIII2,'linewidth',2);
title('Quantized Input Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(xQntIII2)-.1 max(xQntIII2)+.1])
grid

% Calculating Quantized Output
yTemp1III2 = zeros(1,5000);
yQntIII2 = zeros(1,5000); 
yTemp1III2(1) = round((kIII2*xQntIII2(1))*2^(2*DIII2))/2^(2*DIII2);
yQntIII2(1) = round((kIII2*xQntIII2(1))*2^(2*DIII2))/2^(2*DIII2);
yQntIII2(1) = round(yQntIII2(1)*2^(DIII2))/2^(DIII2);
err = zeros(1,5000);

for n = 1:4999    
    err(n) = round((-1*err(n)*2^(2*DIII2)))/(2^(2*DIII2));
    yTemp1III2(n+1) = (kIII2*xQntIII2(n+1))+(aIII2*yQntIII2(n))+err(n);
    yTemp1III2(n+1) = round(yTemp1III2(n+1)*2^(2*DIII2))/2^(2*DIII2);
    yQntIII2(n+1) = round(yTemp1III2(n+1)*2^(DIII2))/2^(DIII2);
    err(n+1) = round((yQntIII2(n+1)-yTemp1III2(n+1))*2^(2*DIII2))/2^(2*DIII2);
end

% Plotting Quantized Output
subplot(223)
plot((0:length(yQntIII2)-1),yQntIII2,'linewidth',2);
title('Quantized Output Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yQntIII2)-.1 max(yQntIII2)+.1])
grid

% Calculating Ideal Output
% As the system transfer function of the noise shaping filter is the same
% as the regular feedback system,
% We can take the ideal system response of the previous filter.
yIdlIII2 = zeros(1,5000);
yIdlIII2(1) = kIII2*xIII2(1);

for n = 1:4999
    yIdlIII2(n+1) = (kIII2*xIII2(n+1))+(aIII2*yIdlIII2(n));
end

% Plotting Ideal Output
subplot(224)
plot((0:length(yIdlIII2)-1),yIdlIII2,'linewidth',2);
title('Ideal Output Signal (D = 8, a = 0.95)')
xlabel('Samples (0-500)')
ylabel('Sample Size')
axis([0 500 -max(yIdlIII2)-.1 max(yIdlIII2)+.1])
grid

yQntErIII2 = yQntIII2-yIdlIII2; % Error between the outputs (Quantization Noise)
VErIII2 = var(yQntErIII2); % Power of the Quantization Noise
VQntSigIII2 = var(yQntIII2); % Power of the Quantized Output Signal
VIdlSigIII2 = var(yIdlIII2); % Power of the Ideal Output SignalSignal

SnRIII2a = 10*log10(VIdlSigIII2/VErIII2);
disp('Ideal Signal to Noise ratio (D = 8, a = 0.95):');disp(SnRIII2a)

SnRIII2b = 10*log10(VQntSigIII2/VErIII2);
disp('Quantization Signal to Noise ratio (D = 8, a = 0.95):');disp(SnRIII2b)

% Theoretical SNR Calculation
theorQNPwrIII2 = 1/(2^(2*DIII2)*12);
theorSnrIII2 = 10*log10(VIdlSigIII2/theorQNPwrIII2);
disp('Theoretical: Ideal Signal to Noise ratio (Noise Shaping Filter) (D = 8, a = 0.95):');disp(theorSnrIII2)

% Plotting the Frequency Responce of the Quantized Output
figure (16)
frqRspQntOPIII2 = 20*log10(fftshift(abs(fft(yQntIII2,4096))));
plot((-29927/2:29927/4096:29927/2 - (29927/4096)),frqRspQntOPIII2, 'b');
title('Frequency Response of the Quantized Output (D = 8, a = 0.95)')
xlabel('Frequency (Hz)')
ylabel('Log Magnitude')
grid

% Explain why the level of the noise is lower at lower frequencies?
% 
% Explanation: The reason behind this is the new noise shaping feedback,
% which pushes the noise on the higher frequencies and acts as a lowpass
% filter. This can be seen from the equation sx(ejw) = [4*((sin(w/2))^2)]*se(ejw)
% At higher frequencies noise power is higher and at lower frequecies noise
% power is lower.
