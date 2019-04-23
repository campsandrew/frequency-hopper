%% General simulation parameters
Rsym = 0.2e6;              % Symbol rate in Hertz
ModulationOrder = 4;       % QPSK alphabet size
Interpolation = 2;         % Interpolation factor
Decimation = 1;            % Decimation factor
Tsym = 1 / Rsym;           % Symbol time in sec
Fs = Rsym * Interpolation; % Sample rate

%% Frame Specifications
BarkerCode = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
BarkerLength = length(BarkerCode);
HeaderLength = BarkerLength * 2;


Message = 'Hello World';
MessageLength = length(Message) + 5;
NumberOfMessage = 100;
PayloadLength = NumberOfMessage * MessageLength * 7;
FrameSize = (HeaderLength + PayloadLength) / log2(ModulationOrder);

%% Tx parameters
RolloffFactor = 0.5;
ScramblerBase = 2;
ScramblerPolynomial = [1 1 1 0 1];
ScramblerInitialConditions = [0 0 0 0];
RaisedCosineFilterSpan = 10;

%% Message generation
msgSet = zeros(100 * MessageLength, 1); 
for msgCnt = 0 : 99
    msgSet(msgCnt * MessageLength + (1 : MessageLength)) = ...
        sprintf('%s %03d\n', Message, msgCnt);
end
integerToBit = comm.IntegerToBit(7, 'OutputDataType', 'double');
MessageBits = integerToBit(msgSet);

% Pluto TX
tx = sdrtx('Pluto', ...
           'CenterFrequency', 915e6, ...
           'BasebandSampleRate', Fs, ...
           'SamplesPerFrame', Interpolation * FrameSize, ...
           'Gain', 0);
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


FrameTime = Interpolation * FrameSize / Fs;
StopTime  = 1000;
currentTime = 0;
while currentTime < StopTime
    data = step(hTx);
    step(tx, data);
    currentTime = currentTime + FrameTime;
end
    
release(hTx);
release(radio);


