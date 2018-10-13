# Quantization-and-Noise-Shaping

Aim of this project is to understand the quantization process and noise shaping technique to increase the SNR which is also called as Sigma-Delta Modulation.

The project is divided into three parts.

In the first part we generate a discrete time signal by sampling a continuous time signal at the sampling frequency fs of 29927 Hz. Then we quantize the signal by rounding the samples to D = 7 bits and find the error between the sampled signal (discrete time signal) and the quantized signal. This error is called as quantization error/noise.
After getting this error signal we find its autocorrelation and plot it. Which turns out to be an impulse proving that the quantization error is nothing but white noise. 
The same tasks are then repeated for D = 12 bits to see how these results change when we increase the number of bits.

In the second part the objective is to find the quantization error at the output of a first order digital filter. For this we first simulate a first order digital filter in MATLAB. Then we again generate a discrete time signal by sampling a continuous time signal at the sampling frequency fs of 29927 Hz. The signal is then quantized with D = 8 bits and applied to the input of the system. We use difference equation method to find the output of the system q(n). We also find the output of the system when the input is not quantized at all, we call this as our ideal output d(n). Then we find the error between d(n) and q(n) which is nothing but the quantization noise. This part also involves the analysis of the SNR and low frequency noise. The same experiment is then repeated for D = 12 bits to see how the results change when we increase the number of bits in our ADC.

Third part involves the use of noise shaping technique in order to increase the SNR (Sigma-Delta Modulation). We repeat all the steps from the the second part but the results are changed we see that now after employing the noise-shaping feed back the noise at lower frequencies is decreased. Then we compared the results of third part where we use noise shaping to increase the SNR with the second part where the system is not using noise shaping.
