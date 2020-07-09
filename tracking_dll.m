settings = initSettings();
channel.PRN = 5;
channel.acquiredFreq = 4.130611212158203e+06;
% channel.codePhase = 16354;
channel.codePhase = 0;
channel.status = 'T';
% [fid, message] = fopen(settings.fileName, 'rb');
[fid, message] = fopen('./myGNSSdata_BB.bin', 'rb');
if (fid > 0)
    % Move the starting point of processing. Can be used to start the
    % signal processing at any point in the data record (e.g. good for long
    % records or for signal processing in blocks).
    fseek(fid, settings.skipNumberOfBytes, 'bof');
else
    fprintf('Cannot open the file!\n');
    return
end
% function [trackResults, channel]= tracking(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute starting 
%                       positions of spreading codes, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: tracking.m,v 1.14.2.32 2007/01/30 09:45:12 dpl Exp $

%% Initialize result structure ============================================
% PID
PID_D = 0;

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);
trackResults.I_E            = zeros(1, settings.msToProcess);
trackResults.I_L            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, settings.msToProcess);
trackResults.Q_P            = zeros(1, settings.msToProcess);
trackResults.Q_L            = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.dllDiscr       = inf(1, settings.msToProcess);
trackResults.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);

% Outputs from LPF
% trackResults.I_P_LPF        = zeros(1, settings.msToProcess * settings.codeLength);
% trackResults.I_E_LPF        = zeros(1, settings.msToProcess * settings.codeLength);
% trackResults.I_L_LPF        = zeros(1, settings.msToProcess * settings.codeLength);
% trackResults.Q_P_LPF        = zeros(1, settings.msToProcess * settings.codeLength);
% trackResults.Q_E_LPF        = zeros(1, settings.msToProcess * settings.codeLength);
% trackResults.Q_L_LPF        = zeros(1, settings.msToProcess * settings.codeLength);

trackResults.I_P_CIC            = zeros(1, settings.msToProcess);

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize filter objects ==========================================
% CIC Decimation Filters
% DO NOT USE: LPF = repmat(cicFilter(settings.samplesPerCode,1,2),1,6);
% Using repmat will create multiple references with the same instance
samplesPerCode = round(settings.samplingFreq/settings.codeFreqBasis*settings.codeLength);
symbolsPerCode = 2;
samplesPerSymbol = 4;
for ii=1:6
    LPF(ii) = cicFilter(1023,1,3);
end

%% Initialize clock sync objects
symbolSync = comm.SymbolSynchronizer(...
    'SamplesPerSymbol',samplesPerSymbol, ...
    'NormalizedLoopBandwidth',0.01, ...
    'DampingFactor',1.0, ...
    'TimingErrorDetector','Gardner (non-data-aided)');

%% Initialize tracking variables ==========================================

codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval
PDIcode = 0.001;

% Calculate filter coefficient values
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);
hwb = waitbar(0,'Tracking...');

%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels
    
    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;
        
        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample) 
        fseek(fid, ...
              settings.skipNumberOfBytes + channel(channelNr).codePhase-1, ...
              'bof');


        % Get a vector with the C/A code sampled 1x/chip
        caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        caCode = [caCode(1023) caCode caCode(1)];

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;

        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;
        preCodeNco   = 0.0;

        %carrier/Costas loop parameters
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;

        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods
            
%% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 50) == 0)
                try
                    plotTracking(1:settings.numberOfChannels, trackResults, settings);
                    waitbar(loopCnt/codePeriods, ...
                            hwb, ...
                            ['Tracking: Ch ', int2str(channelNr), ...
                            ' of ', int2str(settings.numberOfChannels), ...
                            '; PRN#', int2str(channel(channelNr).PRN), ...
                            '; Completed ',int2str(loopCnt), ...
                            ' of ', int2str(codePeriods), ' msec']);                       
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

%% Read next block of data ------------------------------------------------            
            % Find the size of a "block" or code period in whole samples
            
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq;

            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
            
            % Read in the appropriate number of samples to process this
            % interation 
            [rawSignal, samplesRead] = fread(fid, ...
                                             blksize, settings.dataType);
            rawSignal = rawSignal';  %transpose vector
            
            % If did not read in enough samples, then could be out of 
            % data - better exit 
            if (samplesRead ~= blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                fclose(fid);
                return
            end

%% Set up all the code phase tracking information -------------------------
            % Define index into early code vector
            tcode       = (remCodePhase-earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = caCode(tcode2);
            
            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = caCode(tcode2);
            
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = caCode(tcode2);
            
            remCodePhase = (tcode(blksize) + codePhaseStep) - 1023.0;

%% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.samplingFreq;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));
            
            % Finally compute the signal to mix the collected data to bandband
            carrCos = cos(trigarg(1:blksize));
            carrSin = sin(trigarg(1:blksize));

%% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
%             qBasebandSignal = carrCos .* rawSignal;
%             iBasebandSignal = carrSin .* rawSignal;
            qBasebandSignal = rawSignal;
            iBasebandSignal = rawSignal;

            % Now get early, late, and prompt values for each
            I_E = mean(earlyCode  .* iBasebandSignal);
            Q_E = mean(earlyCode  .* qBasebandSignal);
            I_P = mean(promptCode .* iBasebandSignal);
            Q_P = mean(promptCode .* qBasebandSignal);
            I_L = mean(lateCode   .* iBasebandSignal);
            Q_L = mean(lateCode   .* qBasebandSignal);
            
%             if(blksize<=16368)
%                 lessNum = 16368 - blksize;
%                 I_E = LPF(1).CIC_LPF([earlyCode  .* iBasebandSignal,zeros(1,lessNum)]);
%                 Q_E = LPF(2).CIC_LPF([earlyCode  .* qBasebandSignal,zeros(1,lessNum)]);
%                 I_P = LPF(3).CIC_LPF([promptCode .* iBasebandSignal,zeros(1,lessNum)]);
%                 Q_P = LPF(4).CIC_LPF([promptCode .* qBasebandSignal,zeros(1,lessNum)]);
%                 I_L = LPF(5).CIC_LPF([lateCode   .* iBasebandSignal,zeros(1,lessNum)]);
%                 Q_L = LPF(6).CIC_LPF([lateCode   .* qBasebandSignal,zeros(1,lessNum)]);
%             else
%                 fprintf("blksize Error!\r\n");
%                 beep
%             end
            
%             syncOut = symbolSync(I_P_LPF+1j*Q_P_LPF);
            
%             I_E = mean(abs(I_E));
%             I_P = mean(abs(I_P));
%             I_L = mean(abs(I_L));
%             Q_E = mean(abs(Q_E));
%             Q_P = mean(abs(Q_P));
%             Q_L = mean(abs(Q_L));
            
%             I_E = I_E(1);
%             I_P = I_P(1);
%             I_L = I_L(1);
%             Q_E = Q_E(1);
%             Q_P = Q_P(1);
%             Q_L = Q_L(1);
%             
%             if(loopCnt<=2)
%                 I_E = sum(earlyCode  .* iBasebandSignal);
%                 Q_E = sum(earlyCode  .* qBasebandSignal);
%                 I_P = sum(promptCode .* iBasebandSignal);
%                 Q_P = sum(promptCode .* qBasebandSignal);
%                 I_L = sum(lateCode   .* iBasebandSignal);
%                 Q_L = sum(lateCode   .* qBasebandSignal);
%             end
                        
            
%% Generate the six standard LPF output values ---------------------------
            % Now get early, late, and prompt values for each
%             I_E_LPF = earlyCode  .* iBasebandSignal;
%             Q_E_LPF = earlyCode  .* qBasebandSignal;
%             I_P_LPF = promptCode .* iBasebandSignal;
%             Q_P_LPF = promptCode .* qBasebandSignal;
%             I_L_LPF = lateCode   .* iBasebandSignal;
%             Q_L_LPF = lateCode   .* qBasebandSignal;
%% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
%             carrError = atan(Q_P / I_P) / (2.0 * pi);
%             carrError = atan(imag(syncOut(4)) / real(syncOut(4))) / (2.0 * pi);
            carrError = atan(Q_P_LPF(4) / I_P_LPF(4)) / (2.0 * pi);
            
            % Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;
            carrFreq = 0;
            
%             fprintf('D:%f I:%f F:%f\r\n',(tau2carr/tau1carr),(PDIcarr/tau1carr),carrFreqBasis);

            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

%% Find DLL error and update code NCO -------------------------------------
%             codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
%                 (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            codeError = ((I_E * I_E + Q_E * Q_E) - (I_L * I_L + Q_L * Q_L)) / ...
                ((I_E * I_E + Q_E * Q_E) + (I_L * I_L + Q_L * Q_L));
          
            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            codeNco = (codeNco-preCodeNco)*(PID_D) + codeNco;
            preCodeNco = codeNco;
%             fprintf('D:%f I:%f F:%f\r\n',(tau2code/tau1code),(PDIcode/tau1code),settings.codeFreqBasis);
            % Modify code freq based on NCO command
%             codeFreq = settings.codeFreqBasis + codeNco;
            codeFreq = settings.codeFreqBasis;
            
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

%% Record various measures to show in postprocessing ----------------------
            % Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) = ftell(fid);

            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;
            
%             trackResults(channelNr).I_P_CIC(loopCnt) = I_P_CIC;
            
%             N = blksize;
%             bitsRange = ((loopCnt-1)*N+1):loopCnt*N;
%             trackResults(channelNr).I_E_LPF(bitsRange) = I_E_LPF;
%             trackResults(channelNr).I_P_LPF(bitsRange) = I_P_LPF;
%             trackResults(channelNr).I_L_LPF(bitsRange) = I_L_LPF;
%             trackResults(channelNr).Q_E_LPF(bitsRange) = Q_E_LPF;
%             trackResults(channelNr).Q_P_LPF(bitsRange) = Q_P_LPF;
%             trackResults(channelNr).Q_L_LPF(bitsRange) = Q_L_LPF;        
        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;        
        
    end % if a PRN is assigned
end % for channelNr 

% Close the waitbar
close(hwb)
% end
