%% Tx system specs
sampleRateHz = 1e6; % Sample rate
sPS = 2; %samples per symbols
frameSize = 20;
numFrames = 3;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2; %only 4 and 2
filterSymbolSpan = 4;
encoderState=0;
encodedData = [];


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
    mod = comm.BPSKModulator();
    modulatedData = mod.step(data);
    demod = comm.BPSKDemodulator();
else
    mod = comm.QPSKModulator('BitInput', true);
    modulatedData = mod.step(data);
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
freqShift = 0*10^3; 

%% System loop

filteredData = [];%zeros(length(modulatedData)*2,1);
allfilteredTXData = [];
allnoisyData = [];
alloffsetData = [];
allfilteredData = [];
allfilteredDataRef = [];
downsampledRxData = [];
allDownsampledRxData = [];
allDemodulatedData = [];
modtemp = [];

nFD= [];

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
    
    % Time delay signal
%     offsetData = step(varDelay, noisyData, timingOffset); % Variable
    offset= exp(1i*2*pi*freqShift*[1:length(filteredTXData)]); %rx_mystery_waveform.*exp(1i*2*pi*.0195*[1:length(rx_mystery_waveform)]);
    for i = 1:length(filteredTXData)
    offsetData(i) = noisyData(i) *offset(i);
    end
    if k==1
    offsetData= offsetData';
    end
    alloffsetData = [alloffsetData, offsetData];

    % Filter signal
    %filteredData(timeIndex) = step(RxFlt, offsetData);
    nFD = step(RxFlt, noisyData);
    filteredData = [filteredData; nFD]; %#ok<AGROW>
%     filteredDataRef = step(RxFlt, noisyData);
%     allfilteredDataRef = [allfilteredDataRef; filteredDataRef];

    %PLL
    
    %downsample filtered data
    %by Sean Brady
    for i =1:frameSize
        downsampledRxData = [downsampledRxData,filteredData((sPS*filterSymbolSpan)+1+sPS*(i-1))];
    end

    %Demod
    %by Sean Brady
    
    demodulatedData = demod.step(downsampledRxData');
%    allDemodulatedData = [allDemodulatedData,demodulatedData'];
    tdemodulatedData = demod.step(modulatedData);
    % Visualize Error
%     step(cdPre,noisyData);step(cdPost,nFD);pause(0.01); %#ok<*UNRCH>
end
%% Visuals

%testing
x=30
t=1:x;
plot(t,allfilteredTXData(:x)); %filterSymbolSpan*sPS =starting point for rx filter output
hold on;
plot(t, filteredData(5:x+4),'-o','MarkerIndices',1:x);
hold off;
legend('TX','RX');

% cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
%     'Name','Baseband');

% %% PLL setup
% 
% inputDataFull = filteredData;
% InputFrameLen = frameSize;
% 
% %TED
% MaxOutputFrameLen = ceil(InputFrameLen*11/sPS/10);
% alpha = 0.5;
% 
% b = @(Mu) [alpha*Mu^2 - alpha*Mu;... %Interpolator
%     -alpha*Mu^2 - (1-alpha)*Mu + 1;...
%     -alpha*Mu^2 + (1+alpha)*Mu;...
%     alpha*Mu^2 - alpha*Mu];
% 
% LoopFilterState = 0;
% LoopPreviousInput = 0;
% Trigger = false;
% TriggerHistory = false(1, sPS);
% Mu = 0;
% M1Counter = 0;
% InterpFilterState = complex(zeros(3, 1),zeros(3, 1));
% TEDBuffer = complex(zeros(1, sPS), zeros(1, sPS));
% maxOutputSize = ceil(InputFrameLen*11/double(sPS)/10);
% SymbolHolder = complex(zeros(maxOutputSize, 1), zeros(maxOutputSize, 1));
% 
% %Loop Filter
% DetectorGain = 2.7;
% zeta = 1/sqrt(2);
% NormalizedLoopBandwidth = 0.1;
% 
% % Calculate
% Kp = DetectorGain;
% K0 = -1;
% theta = NormalizedLoopBandwidth/sPS/(zeta + 0.25/zeta);
% d  = (1 + 2*zeta*theta + theta^2) * K0 * Kp;
% ProportionalGain = (4*zeta*theta) /d;
% IntegratorGain   = (4*theta*theta)/d;
% 
% %interpolator
% filt = dsp.FIRFilter;
% filt.NumeratorSource='Input port';
% timeCorrect = zeros(size(inputDataFull));
% index = 1;
% buffer = [0 0 0];
% 
% LoopFilter = dsp.IIRFilter( ...
%     'Structure', 'Direct form II transposed', ...
%     'Numerator', [1 0], 'Denominator', [1 -1]);
% NumTrigger  = 0;
% 
% for k=1:frameSize:(numSamples - frameSize)
% 
%     timeIndex = (k:k+frameSize-1).';
%     inputData = inputDataFull(timeIndex); % Grab frame
% 
%     %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Call Synchronizer (sample by sample on frame)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     for i = 1 : InputFrameLen % Process input frame sample-by-sample
% 
%         %% Using Trigger count as output input count
%         NumTrigger = NumTrigger + Trigger;
% 
%         %% Interpolator
%         % Piecewise parabolic interpolator
%         release(filt);
%         % Update coefficents based on new Mu
%         intOut = step(filt,inputData(i),fliplr(b(Mu).'));
% 
%         if Trigger % Interpolation output as symbols
%             SymbolHolder(NumTrigger) = intOut;
%         end
% 
%         %% Timing error detector (TED)
% 
%         % TED calculation occurs on a Trigger
%         if Trigger && all(~TriggerHistory(2:end))
% 
%             % Calculate the midsample from interpolator output
%             t1 = TEDBuffer(end/2 + 1 - rem(SPS,2));
%             t2 = TEDBuffer(end/2 + 1);
%             midSample = (t1+t2)/2;
%             
%             % Data aided method
%             e = real(midSample) * (sign(real(TEDBuffer(1))) - sign(real(intOut))) + ...
%                 imag(midSample) * (sign(imag(TEDBuffer(1))) - sign(imag(intOut)));
%         else
%             e = 0;
%         end
% 
%         % Handle bit stuffing/skipping
%         switch sum([TriggerHistory(2:end), Trigger])
%             case 0
%                 % No Trigger in awhile
%             case 1
%                 % Update buffer (Normal)
%                 TEDBuffer = [TEDBuffer(2:end), intOut];
%             otherwise % > 1
%                 % Stuff condition
%                 TEDBuffer = [TEDBuffer(3:end), 0, intOut];
%         end
% 
%         %% Loop filter
% 
%         v = LoopFilter(e*IntegratorGain);
%         v = e*ProportionalGain + v;
% 
%         %% Interpolation controller
% 
%         W = v + 1/SPS; % W should be small or == SPS when locked
% 
%         TriggerHistory = [TriggerHistory(2:end), Trigger];
% 
%         % Determine if we have an underflow, aka we are at the start edge
%         % of a sample
%         Trigger = (M1Counter < W);
% 
%         % With a new underflow we must update interpolator coefficients
%         if Trigger
%             Mu = M1Counter / W;
%         end
% 
%         % Update counter
%         M1Counter = mod(M1Counter - W, 1);
% 
%     end
% 
%     % Output 1 SPS
%     y = SymbolHolder(1:NumTrigger, 1);
%     timeCorrect(index:index+length(y)-1) = y;
%     index = index + length(y);
% 
% end
% 
% %% Remove zeros
% timeCorrect = timeCorrect(1:index-1);
% 
% %% Visualize
% for k=1:frameSize:(index - frameSize)
%     
%     timeIndex = (k:k+frameSize-1).';
%     step(cdPre,timeCorrect(timeIndex));pause(0.01);
% end
