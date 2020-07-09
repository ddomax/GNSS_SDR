% function genBaseBandData(settings)
%This function generates raw GPS data according 
%to receiver settings for simulation
%
%genData(settings)
%
%   Inputs:
%       settings        - receiver settings.
genSettings.PRN = 5;
genSettings.codeFreq = 1.023e6;
genSettings.remCodePhase = 0.0;
genSettings.IF = settings.IF + 0;
genSettings.bitsPerCode = 1.5;
genSettings.nvRate = (genSettings.codeFreq / settings.codeLength) * genSettings.bitsPerCode;
codePeriods = 1.2 * settings.msToProcess;     % For GPS one C/A code is one ms
signalDuration = codePeriods * 1e-3;
caCode = generateCAcode(genSettings.PRN);
codePhaseStep = genSettings.codeFreq / settings.samplingFreq;
blksize = ceil(settings.codeLength / codePhaseStep);
%% Define index into prompt code vector
tcode       = genSettings.remCodePhase : ...
              codePhaseStep : ...
              ((codePeriods*blksize-1)*codePhaseStep+genSettings.remCodePhase);
tcode2      = ceil(tcode);
tcode3      = mod(tcode2,settings.codeLength) + 1;
promptCode  = caCode(tcode3);
% promptCode  = ones(1,length(promptCode));
%% Generate navigation bits
Ts          = 1 / settings.samplingFreq;
t           = 0:Ts:1.2*(signalDuration-Ts);
nvBitsNum   = 1.3 * signalDuration * genSettings.nvRate;
nvBits      = repmat(generateCAcode(1),1,ceil(nvBitsNum/1023));
bitStep     = genSettings.nvRate/settings.samplingFreq;
bcode       = 0:bitStep:nvBitsNum;
bcode2      = ceil(bcode) + 1;
nvBits2     = nvBits(bcode2);
%% Generate carrier frequency (IF)
LO          = exp(1j*2*pi*genSettings.IF*t);
% baseI       = promptCode .* nvBits2(1:length(promptCode));
baseI       = promptCode;
baseQ       = zeros(1,length(baseI));
baseLen     = length(baseI);
% ifData      = real(LO(1:baseLen)) .* baseI + imag(LO(1:baseLen)) .* baseQ;
ifData      = baseI + baseQ;
ifData      = ifData + 0.00*rand(1,baseLen);
%% Write data file
fid = fopen('./myGNSSdata_BB.bin','wb');
fwrite(fid,int8(ifData*1),settings.dataType);
fclose(fid);
% end