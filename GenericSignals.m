%% by Igor Kavrakov (c)
clear all;close all;clc; matlabrc;
addpath(genpath('CompMet')); %Add the path of the code

% The following script includes an example how to compare two signals. 
% The script is based on the generic signals that are in:
% Kavrakov, I., Kareem, A., and Morgenthal, G. (2020). Comparsion Metric for Time-histories. J. Eng. Mech. 10.1061/(ASCE)EM.1943-7889.0001811
% The signals are described in Eqs. 32-36 in the aforementioned article.

% PLEASE CITE THE ARTICLE IF YOU INTEND TO USE THE METRICS.
% GENERAL USE: If you want to use it for your own signals, you would need
% to define the metric properties (Section 2), based on your own signals
% (Section 1 and 3), and simply call the code as:

% [M]=CompMet(X1,X2,Prop)

% where X1 and X2 are the signals beeing compared (X1 is the reference),
% while Prop is the input structure containing the Metric Properties.
% "M" is the output structure, containing all relevant information on the metrics,
% including all values that are used for plotting. 
% Each call will create a directory with an output name (see below). This will contain
% the M output structure in an .mat file.
% WARNING: Depending on the signal length, M can be a very large file. 

% Copyright (c) Igor Kavrakov, Ahsan Kareem, Guido Morgenthal 2020
%% SECTION 1: General signal Properties (Same for all examples)
dt=0.01; t=dt:dt:100; %Time step/ Time
A=1; A1=1.3;          %Signal amplitude
fc=2;fc1=fc*1.4;fc2=fc*1.8; %Frequencies for 
NSteps=length(t); %Length of signals
SNR=10; %Signal-to-noise ratio
%% SECTION 2: Metric control properties (Here we set the metric properties)
Prop.dt=dt;             %Signal time-step
Prop.TPhase=2/fc;       %Considered significant delay (cf. Eq. 2)

%Metric activision (Which metrics to be activated: 1- Active; 0 - Inactive)
Prop.Metrics.Phase=1;           %Phase & Cross correlation 
Prop.Metrics.Peak=1;            %Peak 
Prop.Metrics.RMS=1;             %RMS 
Prop.Metrics.MagnitudeWarp=1;   %Magnitude Warped 
Prop.Metrics.Magnitude=0;       %Magnitude Unwarped 
Prop.Metrics.Wavelet=1;         %Wavelet 
Prop.Metrics.Stationarity=1;    %Stationarity 
Prop.Metrics.WaveletBicoherence=1; %Bicoherence (nonlinear coupling)
Prop.Metrics.PDF=1;             %Probability Distribution metric 

%METRIC PROPERTIES 
%IMPORTANT: (*) indicates property that is suggested to be modified case-by-case

%Wavelet metric properties
Prop.WavletProperties.nLevel=1200;   %Frequency scales
Prop.WavletProperties.fmax=12;       %(* Maximum frrequency (Strongly dependent on the signals compared)
Prop.WavletProperties.fmin=0.2;       %(* Minimum frequency  (Strongly dependent on the signals compared)
Prop.WavletProperties.beta=3;         %cf. Kijewski & Kareem - beta factor for wavelet (keeping 3 is fine)  
Prop.WavletProperties.f0=6;           %(* Wavelet central frequency (cf. Eq. 10) (Strongly dependent on the signals compared)
Prop.WavletProperties.Padding=0;      %(Recommended 0) Should the signal be padded for the end effects? Generally no
Prop.WavletProperties.MaxDifPlot=.5;  %Plot boundaries (only in case the plots are activated - see below)
Prop.WavletProperties.MaxNormPlot=1.5;%Maximum value for the oclorbar

%Stationarity metric properties
Prop.StationartyProperties.NSurogates=200; %Number of surrogates (Cf. Eq. 16-22)
Prop.StationartyProperties.LocalAnalysis=0; %(Recommended:0) If this is enabled, the local stationarity analysis will be performed regardless of what discriminating statistics says (I.e. Theta in Eq. 14 does not matter).
Prop.StationartyProperties.ConfidenceLevel=0.95; %Confidence level for wavelet surogates Based on Gamma distribution (cf. Eq. 16)
Prop.StationartyProperties.g=2;              %g is exceeding factor for surrogate signals (cf. Eq. 22). 
Prop.StationartyProperties.DiscriminatingStatistic='KLG'; %(Recommended:KLG) Discriminating statistic (KL, LG, KLG) - Kullback–Leibler, Log-Spectral. KLG - mixed.

%Wavelet bispectrum metric properties
Prop.WaveletBicoherenceProperties.Trange=[0 t(end)];     %Time for wavelet analysis (* Dependent on the region of interest)
Prop.WaveletBicoherenceProperties.Rand=1;                %Flag: Rand=1 perform phase randomization. Rand=0 otherwise. (Eq. 26 and 27)
Prop.WaveletBicoherenceProperties.Fact=10*pi;            %Factor for phase randomization (this is r_b for R_b in Eq. 26)
Prop.WaveletBicoherenceProperties.HardTresholding=0;     %(Recommended 0) Hard Tresholding on the bicoherence (Hard value for bicoherence noise)
Prop.WaveletBicoherenceProperties.SoftTresholding=1;     % Building surrogate bispectrum (cf. Eq. 27,28, 29,30,31)
Prop.WaveletBicoherenceProperties.NSurogates=100;        % Number of surrogates for the bispectrum
Prop.WaveletBicoherenceProperties.Eps=0.2;               %Surrogates constant phase interval (this is r_sur near Eq. 26. Ideally should be r_sur=2*pi/r_b)
Prop.WaveletBicoherenceProperties.g=2;                   %Factor of exceedence (cf. Eq. 28)
Prop.WaveletBicoherenceProperties.MaxNormPlot=1.5;       % This is maximum for the plot of the colorbar

%Probability Distribution Metric Properties. Using Kernel method according to Botev (cf. Botev, Z., Grotowski, J., and Kroese, D. (2010). Kernel density estimation via diffusion. Ann.763 Stat., 38, 2916–2957)
Prop.PDFProperties.Kerneldiscretization=2^12; %Discretization of the kernel result (should be power of 2) this is discretization of the integral
Prop.PDFProperties.StandardScore=1;   %Comparison of the Standard Score of the PDFs (i.e. zero mean with normalized standard deviation). Otherwise =0 (full PDF). 
Prop.PDFProperties.MinMax1=0;                  %For the kernel approximation Signal 1- min and max value (WARNING: w.r.t. Standard score, if selected) for approximation (usually this is automatically defined). Used in case of PDF unlimited  values (e.g. noiseless pure sine wave)
Prop.PDFProperties.MinMax2=0;                  %For the kernel approximation Signal 2- min and max value (WARNING: w.r.t. Standard score, if selected) for approximation (usually this is automatically defined)

% PLOTING: properties (Which plots to be activated: 1- Active; 0 - Inactive)
Prop.Plots.General=1;               %Plot time histories & Bar Metrics
Prop.Plots.WarpedSignals=1;         %Plot warped signals
Prop.Plots.Scalogram=1;             %Plot scalograms
Prop.Plots.ScalogramDiff=1;         %Plot difference of scalograms
Prop.Plots.ScalogramNorm=1;         %Plot difference of scalograms (normalised inst freq).
Prop.Plots.MeanPSD=1;               %Plot mean PSD.
Prop.Plots.StatSurrMean=1;          %Plot scalogram filter spectra of surrogates for stationarity analysis
Prop.Plots.StatScalogramFiltered=1; %Plot filtered spectra for stationarity analysis
Prop.Plots.WaveletBicoherence=1;    %Plot wavelet bicoherence
Prop.Plots.WaveletBispectrum=1;     %Plot wavelet bispectrum
Prop.Plots.WaveletBicoherenceFilterd=1;   %Plot wavelet bicoherence filtered with surrogates
Prop.Plots.WaveletBispectrumFilterd=1;   %Plot wavelet bispectrum filtered with surrogates
Prop.Plots.WaveletBispectrumSurogate=1;   %Plot wavelet bispectrum of the surrogates
Prop.Plots.WaveletBispectrumDiff=1;   %Plot wavelet bispectrum of the difference
Prop.Plots.PDFBins=1;              %Plot Optimization of bins (If applicable)
Prop.Plots.PDF=1;                    %Plot PDF
Prop.Plots.FigWidth=16;             %Plot fig width [cm]
Prop.Plots.FigDepth=10;              % Plot fig depth [cm]
%% SECTION 3 CALCULATION 
fprintf(['CompMet: Comparison Metrics for Time-Histories \nIgor Kavrakov, Ahsan Kareem, Guido Morgenthal 2020 (c) \nCite as: Kavrakov, I., Kareem A., and Morgenthal G. (2020). Comparison Metrics for Time-histories: Application to Bridge Aerodynamics. J. Eng. Mech. 10.1061/(ASCE)EM.1943-7889.0001811 \n\n']);
input('!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!\n\nThe following calculation is computationally expensive (due to plotting and the bispectrum metric).\nFor faster computation, turn off the bispectrum metric (Prop.Metrics.WaveletBicoherence=0) and wavelet/stationarity plotting.\n\nPress any key to continue or cancel (CTRL+C).');clc;
%% Test Case 1- Phase shift 
Prop.Name='TestCase1_Phase'; %Name of test case (It will create a directory).

XNoise=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t); %to calculate the STD nosie of the signals
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]); %Generate noise 1
X1=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t+pi/3)+noiseSamp; %Signal 1
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]); %Generate noise signal 2
X2=A.*cos(2.*pi.*fc.*t+pi)+A1.*cos(2.*pi.*fc1.*t+pi/2)+noiseSamp; %Signal 2

CompMet(X1,X2,Prop); %Calculate (X1 - Reference!)
%% Test Case 2- Amplitude scaling
Prop.Name='TestCase2_Amp'; %Name of signal

XNoise=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t); 
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X1=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t+pi/3)+noiseSamp;
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X2=2*A.*cos(2.*pi.*fc.*t)+2*A1.*cos(2.*pi.*fc1.*t+pi/3)+noiseSamp;

CompMet(X1,X2,Prop);
%% Test Case 3- Frequency Modulation (Quadratic chirp);
Prop.Name='TestCase3_Chirp';

K=(fc2-fc1)/t(end); %Frequency modulation factor 
XNoise=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t);
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X1=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t+pi/3)+noiseSamp;
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X2=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*(fc1.*t+K/2.*t.^2)+pi/3)+noiseSamp;

CompMet(X1,X2,Prop);
%% Test Case 3a - Frequency Modulation (Quadratic chirp) - quantification of nonstationary part;
Prop.Name='TestCase3a_ChirpQuant';

K=(fc2-fc1)/t(end);
lamda=0.025;

XNoise=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t);
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X1=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*(fc1.*t+K/2.*t.^2)+pi/3)+noiseSamp;
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X2=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*(fc1.*t+K/2.*t.^2)+pi/3).*exp(-lamda.*t)+noiseSamp;

CompMet(X1,X2,Prop);
%% Test Case 4- Nonlinearity (Quadratic-phase coupling Coupling);
Prop.Name='TestCase4_Coup';

XNoise=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t);
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X1=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t+pi/3)+noiseSamp;
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X2=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t+pi/3)+(A+A1)/2.*cos(2.*pi.*fc.*t).*cos(2.*pi.*fc1.*t+pi/3)+noiseSamp;

CompMet(X1,X2,Prop);
%% Test Case 4a- Nonlinearity+Nonstationarity - quantification of nonstationary part
Prop.Name='TestCase4a_CoupNonstatQuant';

lamda=0.025; %Damping factor for non-stationarity in the nonlinear coupling

XNoise=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t);
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X1=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t+pi/3)+(A+A1)/2.*cos(2.*pi.*fc.*t).*cos(2.*pi.*fc1.*t+pi/3)+noiseSamp;
noiseSamp=random('norm',0,std(XNoise)./SNR,[1 NSteps]);
X2=A.*cos(2.*pi.*fc.*t)+A1.*cos(2.*pi.*fc1.*t+pi/3)+(A+A1)/2.*cos(2.*pi.*fc.*t).*cos(2.*pi.*fc1.*t+pi/3).*exp(-lamda.*t)+noiseSamp;

CompMet(X1,X2,Prop);