% Parameters---------------------------------------------------------------------------------------
% Symbol rate in Hertz
Rsym = 0.2e6;
% QPSK alphabet size
ModulationOrder = 4; 
% Interpolation factor
Interpolation = 2;
% Decimation factor
Decimation = 1;
% Symbol time in sec
Tsym = 1/Rsym;  
% Sample rate
Fs = Rsym * Interpolation; 
channels = [2.4e9, 2.401e9, 2.402e9,2.403e9,2.404e9,2.405e9,2.406e9,2.407e9,2.408e9];
%Reciever
%---------------------------------------------------------------------------------------------------
% Frame Specifications
% Bipolar Barker Code
BarkerCode = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];    
BarkerLength = length(BarkerCode);
HeaderLength = BarkerLength * 2;                   
Message = 'ECE 531 Project Message';
MessageLength = length(Message) + 5;               
NumberOfMessage = 100;  
% 7 bits per characters
PayloadLength = NumberOfMessage * MessageLength * 7;
% Frame size in symbols
FrameSize = (HeaderLength + PayloadLength)/ log2(ModulationOrder);
FrameTime = Tsym*FrameSize;
%Rx parameters
RolloffFactor = 0.5;                     
ScramblerBase = 2;
ScramblerPolynomial = [1 1 1 0 1];
ScramblerInitialConditions = [0 0 0 0];
RaisedCosineFilterSpan = 10;                  
DesiredPower = 2;            
AveragingLength = 50;           
MaxPowerGain = 60;           
MaximumFrequencyOffset = 6e3;
K = 1;
A = 1/sqrt(2);
PhaseRecoveryLoopBandwidth = 0.01;         
PhaseRecoveryDampingFactor = 1;            
TimingRecoveryLoopBandwidth = 0.01;         
TimingRecoveryDampingFactor = 1;            
TimingErrorDetectorGain = 2.7*2*K*A^2+2.7*2*K*A^2;
PreambleDetectorThreshold = 0.8;
% Message generation and BER calculation parameters
msgSet = zeros(100 * MessageLength, 1);
for msgCnt = 0 : 99
    msgSet(msgCnt * MessageLength + (1 : MessageLength)) = ...
        sprintf('%s %03d\n', Message, msgCnt);
end
integerToBit = comm.IntegerToBit(7, 'OutputDataType', 'double');
MessageBits = integerToBit(msgSet);
% For BER calculation masks
BerMask = zeros(NumberOfMessage * length(Message) * 7, 1);
for i = 1 : NumberOfMessage
    BerMask( (i-1) * length(Message) * 7 + ( 1: length(Message) * 7) ) = ...
        (i-1) * MessageLength * 7 + (1: length(Message) * 7);
end
% Pluto receiver parameters
PlutoCenterFrequency = channels(1);
PlutoGain = 30;
PlutoFrontEndSampleRate = Fs;
PlutoFrameLength = Interpolation * FrameSize;
% Experiment parameters
PlutoFrameTime = PlutoFrameLength / PlutoFrontEndSampleRate;
StopTime = 1000;
hRx  = QPSKReceiver(...
    'ModulationOrder', ModulationOrder, ...
    'SampleRate', Fs, ...
    'DecimationFactor', Decimation, ...
    'FrameSize', FrameSize, ...
    'HeaderLength', HeaderLength, ...
    'NumberOfMessage', NumberOfMessage, ...
    'PayloadLength', PayloadLength, ...
    'DesiredPower', DesiredPower, ...
    'AveragingLength', AveragingLength, ...
    'MaxPowerGain', MaxPowerGain, ...
    'RolloffFactor', RolloffFactor, ...
    'RaisedCosineFilterSpan', RaisedCosineFilterSpan, ...
    'InputSamplesPerSymbol', Interpolation, ...
    'MaximumFrequencyOffset', MaximumFrequencyOffset, ...
    'PostFilterOversampling', Interpolation/Decimation, ...
    'PhaseRecoveryLoopBandwidth', PhaseRecoveryLoopBandwidth, ...
    'PhaseRecoveryDampingFactor', PhaseRecoveryDampingFactor, ...
    'TimingRecoveryDampingFactor', TimingRecoveryDampingFactor, ...
    'TimingRecoveryLoopBandwidth', TimingRecoveryLoopBandwidth, ...
    'TimingErrorDetectorGain', TimingErrorDetectorGain, ...
    'PreambleDetectorThreshold', PreambleDetectorThreshold, ...
    'DescramblerBase', ScramblerBase, ...
    'DescramblerPolynomial', ScramblerPolynomial, ...
    'DescramblerInitialConditions', ScramblerInitialConditions,...
    'BerMask', BerMask,...
    'PrintOption',true);
%Transmitter
%-------------------------------------------------------------------------------------------------
hTx = QPSKTransmitter(...
    'UpsamplingFactor',             Interpolation, ...
    'RolloffFactor',                RolloffFactor, ...
    'RaisedCosineFilterSpan',       RaisedCosineFilterSpan, ...
    'MessageBits',                  MessageBits, ...
    'MessageLength',                MessageLength, ...
    'NumberOfMessage',              NumberOfMessage, ...
    'ScramblerBase',                ScramblerBase, ...
    'ScramblerPolynomial',          ScramblerPolynomial, ...
    'ScramblerInitialConditions',   ScramblerInitialConditions);
% Frame Specifications
BarkerCode = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
BarkerLength = length(BarkerCode);
HeaderLength = BarkerLength * 2;
Message = 'Freq: 1';
MessageLength = length(Message) + 5;
NumberOfMessage = 100;
PayloadLength = NumberOfMessage * MessageLength * 7;
FrameSize = (HeaderLength + PayloadLength) / log2(ModulationOrder);
% Tx parameters
RolloffFactor = 0.5;
ScramblerBase = 2;
ScramblerPolynomial = [1 1 1 0 1];
ScramblerInitialConditions = [0 0 0 0];
RaisedCosineFilterSpan = 10;
% Pluto transmit and recieve objects
txp = sdrtx('Pluto');
txp.CenterFrequency       = channels(1);
txp.BasebandSampleRate    = PlutoFrontEndSampleRate;
txp.SamplesPerFrame       = PlutoFrameLength;
txp.GainSource            = 'Manual';
txp.Gain                  = 0;
txp.OutputDataType        = 'double';
rxp = sdrrx('Pluto');
rxp.CenterFrequency       = channels(1);
rxp.BasebandSampleRate    = PlutoFrontEndSampleRate;
rxp.SamplesPerFrame       = PlutoFrameLength;
rxp.GainSource            = 'Manual';
rxp.Gain                  = 30;
rxp.OutputDataType        = 'double';
% Initialize variables
currentTime = 0;
BER = [];
rcvdSignal = complex(zeros(PlutoFrameLength,1));
curr_freq = 1;
BER_T = 1;
while currentTime < StopTime
    %Receive signal from the radio
    rcvdSignal = rxp();
    %Decode the received message
    [~, ~, ~, BER] = hRx(rcvdSignal);
    BER_C = BER(1);
    %Decide on changing frequency
    if BER_C > BER_T
        curr_freq = randi([1,11]);
        txp.CenterFrequency = channels(curr_freq);
        rxp.CenterFrequency = channels(curr_freq);
    end
    %Send frequency over QPSK
    tx(hTx(currFreq));
    %Update simulation time
    currentTime=currentTime+(radio.SamplesPerFrame / radio.BasebandSampleRate);
end
release(htx);
release(rx);
release(txp);
release(rxp);
