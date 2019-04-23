%% Variable Initializations
visuals = true;                                         % Debug flag
loFreq = 2.4e9;                                         % LO Center Frequency
sampleRate = 1e6;
samplesPerSymbol = 2;
filterSymbolSpan = 4;
modOrder = 2;
frameSize = 2^10;
numFrames = 100;
randData = randi([0 1], [numFrames, frameSize]);        % Binary data

% Impairments Variables
snr = 20;                                               % Noise (SNR)
timingOffset = samplesPerSymbol * 0.01;                 % Timing Offset
freqOffset = sampleRate * 0.1;                          % Frequency Offset
phaseOffset = 0;                                        % Phase Offset

% Object Variables
varDelay = dsp.VariableFractionalDelay();
dbpskMod = comm.DBPSKModulator();                       % Modulator
dbpskDemod = comm.DBPSKDemodulator();                   % Demodulator
txFlt = comm.RaisedCosineTransmitFilter(...             % TX Filter
    'OutputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan);
rxFlt = comm.RaisedCosineReceiveFilter(...              % RX Filter
    'InputSamplesPerSymbol', samplesPerSymbol,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', 1);
chan = comm.AWGNChannel('NoiseMethod', ...              % AWGN Channel
    'Signal to noise ratio (SNR)', ...
    'SNR', snr, ...
    'SignalPower', 1, ...
    'RandomStream', 'mt19937ar with seed');
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
pfOffset = comm.PhaseFrequencyOffset('FrequencyOffset', freqOffset, ...
    'PhaseOffset', phaseOffset, ...
    'SampleRate', sampleRate);

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

%% Start Simulation
for frame = 1:numFrames
    dataFrame = randData(frame, :).';                       % Data
    txModData = dbpskMod(dataFrame);                        % Modulate Data
    txFiltData = txFlt(txModData);                          % TX Filter Data    
    
    %TODO: Check for receive message based on
    %      receiver message bit error
%     
    % Add imparments to data
    txData = varDelay(txFiltData, frame * timingOffset);    % Timing Offset Data
    txData = pfOffset(txData);                              % Phase Frequency Data
    txData = chan(txData);                                  % Noisy Data
%     disp("Center Frequency: " + loFreq);
    
    % Create random channel usage
    %TODO: Add channel usage
    
    % Matched Filter and Corrections
    rxFiltData = rxFlt(txData);                    % Receiver Filter
    rxData = coarseSync(rxFiltData);                   % Coarse Frequency Correction
    rxData = symbolSync(rxData);               % Timing Correction
    tSize = size(rxData);                          % Correct frame sizing after timing correction
    if tSize(1) < frameSize
        stuff = frameSize - tSize(1);
        rxData = [rxData; zeros(stuff, 1)];
    end
    if tSize(1) > frameSize
       remove = tSize(1) - frameSize;
       rxData = rxData(1:end - remove);
    end
    rxData = fineSync(rxData);                     % Fine Frequency Correction 
    %TODO: Frame syncronization
    rxDemodData = dbpskDemod(rxData);          % Demodulate Data
    
    
    %TODO: Calculate bit error
    %      Check if bit error is too high
    %      Create random frequency to change to
    %      Transmit frequency change back    
%     e = biterr(dataFrame, rxDemodData);
%     disp(e);
    
    % Frame by frame graphs
    if visuals
        cdRef(txFiltData);              % After TX Filter
        cdPre(txData);                  % Sent Data
        cdPost(rxData);                 % Data After Corrections
        pause(0.1);
    end
end