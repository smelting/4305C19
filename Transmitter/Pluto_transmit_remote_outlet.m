bitsPerSecond = floor(1.0/0.00018195);  %from baseband timing

upsampleRate = 150;
sampleRateHz = bitsPerSecond*upsampleRate*4.0; % Sample rate
recordTime = .05; %10 secsonds of recording
samplesPerFrame = sampleRateHz*recordTime;
centerFreq = 433.920e+06;

data_chan3_on =   [0;1;0;1;0;0;0;0;0;1;0;1;0;1;1;1;0;0;0;0;0;0;1;1;0;];
data_chan3_off =  [0;1;0;1;0;0;0;0;0;1;0;1;0;1;1;1;0;0;0;0;1;1;0;0;0;];
data_chan2_on =   [0;1;0;1;0;0;0;0;0;1;0;1;0;1;0;1;1;1;0;0;0;0;1;1;0;];
data_chan2_off =  [0;1;0;1;0;0;0;0;0;1;0;1;0;1;0;1;1;1;0;0;1;1;0;0;0;];
data_chan1_on =   [0;1;0;1;0;0;0;0;0;1;0;1;0;1;0;1;0;0;1;1;0;0;1;1;0;];
data_chan1_off =  [0;1;0;1;0;0;0;0;0;1;0;1;0;1;0;1;0;0;1;1;1;1;0;0;0;];

%%Setup Visuals
sa = dsp.SpectrumAnalyzer('SampleRate',sampleRateHz,'ShowLegend',true);
ts = dsp.TimeScope('SampleRate',sampleRateHz,'TimeSpan',recordTime);

cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');

transmit_chan3_on = Encoded_Remote(data_chan3_on,upsampleRate,2000);
transmit_chan3_off = Encoded_Remote(data_chan3_off,upsampleRate,2000);
transmit_chan2_on = Encoded_Remote(data_chan2_on,upsampleRate,2000);
transmit_chan2_off = Encoded_Remote(data_chan2_off,upsampleRate,2000);
transmit_chan1_on = Encoded_Remote(data_chan1_on,upsampleRate,2000);
transmit_chan1_off = Encoded_Remote(data_chan1_off,upsampleRate,2000);

%transmit = complex(upsampled,-1*upsampled);
%%Setup Radio
rx =sdrrx('Pluto','OutputDataType','double','GainSource','AGC Fast Attack','SamplesPerFrame',samplesPerFrame,'CenterFrequency',centerFreq,'BasebandSampleRate',sampleRateHz)
tx = sdrtx('Pluto','Gain',0,'CenterFrequency',centerFreq,'SamplesPerFrame',samplesPerFrame,'BasebandSampleRate',sampleRateHz)


%This seems to make the tx a big less ugly 
TxFlt = comm.RaisedCosineTransmitFilter(... 
    'OutputSamplesPerSymbol', 4,...
    'FilterSpanInSymbols', upsampleRate*8);

output_chan3_on = step(TxFlt,transmit_chan3_on);
output_chan3_off = step(TxFlt,transmit_chan3_off);
output_chan2_on = step(TxFlt,transmit_chan2_on);
output_chan2_off = step(TxFlt,transmit_chan2_off);
output_chan1_on = step(TxFlt,transmit_chan1_on);
output_chan1_off = step(TxFlt,transmit_chan1_off);

tx.transmitRepeat(output_chan3_off);


ts.release();
while 1
%    [rxData,validData,overflow] = rx();
%    ts(rxData);
%    sa(rxData);
%    step(cdPre,rxData);
 disp("Turning on channel 3");
 tx.release();
 tx.transmitRepeat(output_chan1_on);
 pause(3);
 tx.release();
 tx.transmitRepeat(output_chan3_on);
 pause(3);
 disp("Turning off channel 3");
 tx.release();
 tx.transmitRepeat(output_chan1_off);
 pause(3);
 tx.release();
 tx.transmitRepeat(output_chan3_off);
 pause(3);

end