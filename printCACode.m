clc;
%--- Include folders with functions ---------------------------------------
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions

%
caCode = generateCAcode(22); % PRN=5 for myGNSSdata_xxxx.bin PRN=22 for gioveA&B

fprintf('[');
for ii=1:(length(caCode)-1)
    fprintf('%d,',caCode(ii));
end
fprintf('%d]',caCode(length(caCode)));