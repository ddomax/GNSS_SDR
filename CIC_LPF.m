function y = CIC_LPF(x)
%CIC_LPF Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 9.2 and the DSP System Toolbox 9.4.
% Generated on: 01-Jul-2020 15:03:14

%#codegen

% To generate C/C++ code from this function use the codegen command.
% Type 'help codegen' for more information.
x = x';

persistent Hd;

if isempty(Hd)
    
    decf    = 16368;  % Decimation Factor
    diffd   = 2;      % Differential Delay
    numsecs = 6;      % Number of Sections
    
    Hd = dsp.CICDecimator( ...
        'DecimationFactor', decf, ...
        'DifferentialDelay', diffd, ...
        'NumSections', numsecs);
end

s = fi(x,1,32,28,'RoundingMethod','Round','OverflowAction','Saturate');
y_fix = step(Hd,s);
% y = (y_fix.data)./16367.5004882813; %diifd1,secs2
y = (y_fix.data)./7.51896487170364e+22; %diffd2,secs6
% y = y_fix.data;


% [EOF]
