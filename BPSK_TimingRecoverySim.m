close all;
M = 2;         % Modulation order for BPSK
nSym = 20000;  % Number of symbols in a packet
sps = 4;       % Samples per symbol
timingErr = 2; % Samples of timing error
snr = 15;      % Signal-to-noise ratio (dB)

txfilter = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol',sps);
rxfilter = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol',sps,'DecimationFactor',1);

symbolSync = comm.SymbolSynchronizer(...
    'SamplesPerSymbol',sps, ...
    'NormalizedLoopBandwidth',0.01, ...
    'DampingFactor',1.0, ...
    'TimingErrorDetector','Gardner (non-data-aided)');

data = randi([0 M-1],nSym,1);
modSig = pskmod(data,M);

fixedDelay = dsp.Delay(timingErr);
fixedDelaySym = ceil(fixedDelay.Length/sps); % Round fixed delay to nearest integer in symbols

txSig = txfilter(modSig);
delayedSig = fixedDelay(txSig);

rxSig = awgn(delayedSig,snr,'measured');

rxSample = rxfilter(rxSig);
scatterplot(rxSample(10000:end),2)

rxSync = symbolSync(rxSample);
scatterplot(rxSync(10000:end),2)

recData = pskdemod(rxSync,M);

sysDelay = dsp.Delay(fixedDelaySym + txfilter.FilterSpanInSymbols/2 + ...
    rxfilter.FilterSpanInSymbols/2);

[numErr1,ber1] = biterr(sysDelay(data),recData(1:length(data)))