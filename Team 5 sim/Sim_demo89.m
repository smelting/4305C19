%% Tx system specs
sampleRateHz = 1e6; % Sample rate
sPS = 4; %samples per symbols
frameSize = 73;
numFrames = 1;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2; %only 4 and 2
filterSymbolSpan = 4;
encoderState=0;
encodedData = [];

%% Packet Structure

%9 bit preamble
%32b ID
%8b temp
%8b pressure
%8b flags= error message
%8b CRC = check for data corruption, need to find in real life
%73b frameSize


%% Data setup
data = randi([0 1], numSamples, 1);

%% Encoder

if modulationOrder== 4
    tempData = [];
   for i= 1:2:length(data)
       tempData = [data(i);data(i+1)];
       encodedData = [encodedData, tempData];
   end
end
encodedData = encodedData';
%% Modulation
if modulationOrder == 2
    modl = comm.BPSKModulator();
    modulatedData = modl.step(data);
    demod = comm.BPSKDemodulator();
else
    modl = comm.QPSKModulator('BitInput', true);
    modulatedData = modl.step(data);
    demod = comm.QPSKDemodulator();
end
%% Matched Filter Setup
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', sPS,...
    'FilterSpanInSymbols', filterSymbolSpan);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', sPS,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', 1) ;% Set to filterUpsample/2 when introducing timing estimation

%% Frequency Distortion and AWGN setup
noiseStrength = 0;
freqShift = 2*10^3; 

%% PLL setup

%TED
MaxOutputFrameLen = ceil(frameSize*11/sPS/10);
alpha = 0.5;

b = @(Mu) [alpha*Mu^2 - alpha*Mu;... %Interpolator
    -alpha*Mu^2 - (1-alpha)*Mu + 1;...
    -alpha*Mu^2 + (1+alpha)*Mu;...
    alpha*Mu^2 - alpha*Mu];

LoopFilterState = 0;
LoopPreviousInput = 0;
Trigger = false;
TriggerHistory = false(1, sPS);
Mu = 0;
M1Counter = 0;
InterpFilterState = complex(zeros(3, 1),zeros(3, 1));
TEDBuffer = complex(zeros(1, sPS), zeros(1, sPS));
maxOutputSize = ceil(frameSize*11/double(sPS)/10);
SymbolHolder = complex(zeros(maxOutputSize, 1), zeros(maxOutputSize, 1));

%Loop Filter
DetectorGain = 2.7;
zeta = 1/sqrt(2);
NormalizedLoopBandwidth = 0.1;

% Calculate
Kp = DetectorGain;
K0 = -1;
theta = NormalizedLoopBandwidth/sPS/(zeta + 0.25/zeta);
d  = (1 + 2*zeta*theta + theta^2) * K0 * Kp;
ProportionalGain = (4*zeta*theta) /d;
IntegratorGain   = (4*theta*theta)/d;

%interpolator
filt = dsp.FIRFilter;
filt.NumeratorSource='Input port';
%timeCorrect = zeros(size(filteredData));
index = 1;
buffer = [0 0 0];

LoopFilter = dsp.IIRFilter( ...
    'Structure', 'Direct form II transposed', ...
    'Numerator', [1 0], 'Denominator', [1 -1]);



%% System loop

filteredData = [];%zeros(length(modulatedData)*2,1);
allfilteredTXData = [];
allnoisyData = [];
alloffsetData = [];
allfilteredData = [];
downsampledRxData = [];
allDownsampledRxData = [];
allDemodulatedData = [];
modtemp = [];

nFD= [];
allintOut= [];

NumTrigger = 0;

pllData = [];

pll = comm.SymbolSynchronizer(...
    'SamplesPerSymbol', sPS, ...
    'NormalizedLoopBandwidth', 0.01, ...
    'DampingFactor', 1.0, ...
    'DetectorGain', 20000, ...
    'TimingErrorDetector', 'Zero-Crossing (decision-directed)');

for k=1:frameSize:(numSamples)
    
    timeIndex = (k:k+frameSize-1).';
    
    % Filter signal
    for fltfix= 1:sPS
    modtemp = [modulatedData(timeIndex);0;0;0;0]; %4
    end
    filteredTXData = step(TxFlt, modtemp);
    for fltfix= 1:(sPS*filterSymbolSpan)
    filteredTXData =[filteredTXData; 0;]; %10
    end
    allfilteredTXData = [allfilteredTXData; filteredTXData];
    
    % Pass through channel
    noiseChan = noiseStrength*randn([1 length(filteredTXData)]); %%AWGN
    noisyData = filteredTXData + noiseChan';
    allnoisyData = [allnoisyData; noisyData];
    
    % Filter signal
    %filteredData(timeIndex) = step(RxFlt, offsetData);
    nFD = step(RxFlt, filteredTXData);
    filteredData = [filteredData; nFD]; %#ok<AGROW>
    %downsample filtered data
    %by Sean Brady
%     for i =1:frameSize
%         downsampledRxData = [downsampledRxData,filteredData((sPS*filterSymbolSpan)+1+sPS*(i-1))];
%     end

    %Demod
    %by Sean Brady
   % demodulatedData = demod.step(downsampledRxData');
%    allDemodulatedData = [allDemodulatedData,demodulatedData'];
    % Visualize Error
%     step(cdPre,noisyData);step(cdPost,nFD);pause(0.01); %#ok<*UNRCH>
end

%testing
t= 1:73;
plot(t,filteredTXData(1:73));
hold on;
plot(t,filteredData(9:81))
hold off;

