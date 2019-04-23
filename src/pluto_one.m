%% Variable Initializations
visuals = true;                                         % Debug flag
loFreq = 2.4e9;                                         % LO Center Frequency
sampleRate = 1e6;
samplesPerSymbol = 2;
filterSymbolSpan = 4;
modOrder = 1;
frameSize = 2^10;

% Object Variables
dbpskMod = comm.DBPSKModulator();                       % Modulator
dbpskDemod = comm.DBPSKDemodulator();                   % Demodulator
rx = sdrrx('Pluto',...                                  % Pluto Receiver
    'CenterFrequency', loFreq, ...
    'SamplesPerFrame', frameSize, ...    
    'OutputDataType', 'double', ...
    'BasebandSampleRate', sampleRate);
tx = sdrtx('Pluto', ...                                 % Pluto Transmitter
    'CenterFrequency', loFreq, ...
    'BasebandSampleRate', sampleRate);
txFlt = comm.RaisedCosineTransmitFilter(...             % TX Filter
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan);
symbolSync = comm.SymbolSynchronizer(...
    'SamplesPerSymbol', samplesPerSymbol, ...
    'NormalizedLoopBandwidth', 0.01, ...
    'DampingFactor', 1.0, ...
    'TimingErrorDetector','Zero-Crossing (decision-directed)');
coarseSync = comm.CoarseFrequencyCompensator('Modulation','BPSK', ...
    'FrequencyResolution', 1, ...
    'SampleRate', sampleRate * samplesPerSymbol);
fineSync = comm.CarrierSynchronizer('DampingFactor', 0.7, ...
    'NormalizedLoopBandwidth', 0.005, ...
    'SamplesPerSymbol', samplesPerSymbol, ...
    'Modulation','BPSK');

% Visual Objects
% sa = dsp.SpectrumAnalyzer('SampleRate', sampleRate, 'ShowLegend', true);
cdRef = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Reference');
cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband with Timing Offset');
cdPost = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Corrected Timing Offset');
cdRef.Position(1) = 50;
cdPre.Position(1) = cdRef.Position(1) + cdRef.Position(3) + 10;
cdPost.Position(1) = cdPre.Position(1) + cdPre.Position(3) + 10;

%% Transmit Repeat Useful Data
% tx.transmitRepeat();

%% Receive Data
while true
    rxData = rxFlt(rx());                          % Receive data
    
    % Apply corrections
    rxData = symbolSync(rxFiltData);               % Timing Correction
    tSize = size(rxData);                          % Correct frame sizing after timing
    if tSize(1) < frameSize
        stuff = frameSize - tSize(1);
        rxData = [rxData; zeros(stuff, 1)];
    end
    if tSize(1) > frameSize
       remove = tSize(1) - frameSize;
       rxData = rxData(1:end - remove);
    end
    rxData = coarseSync(rxData);                   % Coarse Frequency Correction
    rxData = fineSync(rxData);                     % Fine Frequency Correction 
    %TODO: Frame syncronization
    rxDemodData = dbpskDemod(rxFiltData);              % Demodulate Data
    
    % Frame by frame graphs
    if visuals
%         cdRef(txFiltData);              % After TX Filter
%         cdPre(txData);                  % Sent Data
        cdPost(rxData);                 % Data After Corrections
%         pause(0.1);
    end 
end