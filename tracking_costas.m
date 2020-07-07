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
%% Clean up the environment first =========================================
clear; close all; clc;

format ('compact');
format ('long', 'g');

%--- Include folders with functions ---------------------------------------
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions

%---- Init settings
settings = initSettings();
settings.fileName           = ...
   '.\myGNSS_NO_PN.bin';
%% Initialize result structure ============================================
% PID
PID_D = 0;

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, settings.msToProcess);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, settings.msToProcess);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_P            = zeros(1, settings.msToProcess);

% Loop discriminators
trackResults.pllDiscr       = inf(1, settings.msToProcess);
trackResults.pllDiscrFilt   = inf(1, settings.msToProcess);

%% Initialize filter objects ==========================================
% CIC Decimation Filters
% DO NOT USE: LPF = repmat(cicFilter(settings.samplesPerCode,1,2),1,6);
% Using repmat will create multiple references with the same instance
samplesPerCode = round(settings.samplingFreq/settings.codeFreqBasis*settings.codeLength);
symbolsPerCode = 2;
samplesPerSymbol = 4;
DecimationRatio = 8;
for ii=1:2
    LPF(ii) = cicFilter(DecimationRatio,1,3);
end

%% Initialize clock sync objects
symbolSync = comm.SymbolSynchronizer(...
    'SamplesPerSymbol',samplesPerSymbol, ...
    'NormalizedLoopBandwidth',0.01, ...
    'DampingFactor',1.0, ...
    'TimingErrorDetector','Gardner (non-data-aided)');

%% Initialize tracking variables ==========================================
blksize = DecimationRatio;
codePeriods = 1e-3*settings.msToProcess*settings.samplingFreq/blksize;     % For GPS one C/A code is one ms

%--- PLL variables --------------------------------------------------------
% Summation interval
% PDIcarr = 1/(settings.samplingFreq/blksize);
PDIcarr = 0.001;

% Calculate filter coefficient values
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);
hwb = waitbar(0,'Tracking...');

%% Start processing channels ==============================================
% Move the starting point of processing. Can be used to start the
% signal processing at any point in the data record (e.g. for long
% records). In addition skip through that data file to start at the
% appropriate sample (corresponding to code phase). Assumes sample
% type is schar (or 1 byte per sample) 
[fid, message] = fopen(settings.fileName, 'rb');
fseek(fid, ...
      settings.skipNumberOfBytes, ...
      'bof');

%--- Perform various initializations ------------------------------

% define initial code frequency basis of NCO
codeFreq      = settings.codeFreqBasis;

% define carrier frequency which is used over whole tracking period
% carrFreq      = settings.IF-100;
carrFreq        = 0;
carrFreqBasis = settings.IF-100;
% define residual carrier phase
remCarrPhase  = 0.0;

%carrier/Costas loop parameters
oldCarrNco   = 0.0;
oldCarrError = 0.0;

%=== Process the number of specified code periods =================
for loopCnt =  1:codePeriods

%% GUI update -------------------------------------------------------------
    % The GUI is updated every 50ms. This way Matlab GUI is still
    % responsive enough. At the same time Matlab is not occupied
    % all the time with GUI task.
    if (rem(loopCnt, 5000) == 0)
        try
            plotTracking_costas(trackResults, settings);
            waitbar(loopCnt/codePeriods, ...
                    hwb, ...
                    ['Completed ',int2str(loopCnt), ...
                    ' of ', int2str(codePeriods), ' blocks']);                       
        catch
            % The progress bar was closed. It is used as a signal
            % to stop, "cancel" processing. Exit.
            disp('Progress bar closed, exiting...');
            return
        end
    end

%% Read next block of data ------------------------------------------------            
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
    qBasebandSignal = carrCos .* rawSignal;
    iBasebandSignal = carrSin .* rawSignal;

    % Now get early, late, and prompt values for each
    I_P = sum(iBasebandSignal);
    Q_P = sum(qBasebandSignal);

        I_P_LPF = LPF(1).CIC_LPF(iBasebandSignal);
        Q_P_LPF = LPF(2).CIC_LPF(qBasebandSignal);


%             syncOut = symbolSync(I_P_LPF+1j*Q_P_LPF);

%% Find PLL error and update carrier NCO ----------------------------------

    % Implement carrier loop discriminator (phase detector)
%     carrError = atan(Q_P / I_P) / (2.0 * pi);
%     carrError = atan(imag(syncOut(4)) / real(syncOut(4))) / (2.0 * pi);
%     if(I_P_LPF==0)
%         I_P_LPF = 0.1;
%     end
%     carrError = atan(Q_P_LPF / I_P_LPF) / (2.0 * pi);
    carrError = sign(I_P_LPF)*angle(abs(I_P_LPF)+1j*Q_P_LPF) / (2.0 * pi);

    % Implement carrier loop filter and generate NCO command
    carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
        (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
    oldCarrNco   = carrNco;
    oldCarrError = carrError;

    % Modify carrier freq based on NCO command
    carrFreq = carrFreqBasis + carrNco;

%     fprintf('D:%f I:%f F:%f\r\n',(tau2carr/tau1carr),(PDIcarr/tau1carr),carrFreqBasis);

    trackResults.carrFreq(loopCnt) = carrFreq;

%% Record various measures to show in postprocessing ----------------------
    % Record sample number (based on 8bit samples)
    trackResults.absoluteSample(loopCnt) = ftell(fid);

    trackResults.pllDiscr(loopCnt)       = carrError;
    trackResults.pllDiscrFilt(loopCnt)   = carrNco;

    trackResults.I_P(loopCnt) = I_P;
    trackResults.Q_P(loopCnt) = Q_P;

end % for loopCnt    
        
plotTracking_costas(trackResults, settings)

% Close the waitbar
close(hwb)
