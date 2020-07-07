function plotTracking_costas(trackResults, settings)
%This function plots the tracking results for the given channel list.
%
%plotTracking(channelList, trackResults, settings)
%
%   Inputs:
%       channelList     - list of channels to be plotted.
%       trackResults    - tracking results from the tracking function.
%       settings        - receiver settings.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
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
%$Id: plotTracking.m,v 1.5.2.23 2006/08/14 14:45:14 dpl Exp $

%=== For all listed channels ==============================================

%% Select (or create) and clear the figure ================================
    % The number 200 is added just for more convenient handling of the open
    % figure windows, when many figures are closed and reopened.
    % Figures drawn or opened by the user, will not be "overwritten" by
    % this function.

    figure(200);
    clf(200);

%% Draw axes ==============================================================
        % Row 1
        handles(1, 1) = subplot(3, 3, 1);
        handles(1, 2) = subplot(3, 3, [2 3]);
        % Row 2
        handles(2, 1) = subplot(3, 3, 4);
        handles(2, 2) = subplot(3, 3, [5 6]);
        % Row 3
        handles(3, 1) = subplot(3, 3, 7);
        handles(3, 2) = subplot(3, 3, 8);
        handles(3, 3) = subplot(3, 3, 9);

%% Plot all figures =======================================================

        timeAxisInSeconds = (1:settings.msToProcess)/1000;

        %----- Discrete-Time Scatter Plot ---------------------------------
        plot(handles(1, 1), trackResults.I_P,...
                            trackResults.Q_P, ...
                            '.');

        grid  (handles(1, 1));
        axis  (handles(1, 1), 'equal');
        title (handles(1, 1), 'Discrete-Time Scatter Plot');
        xlabel(handles(1, 1), 'I prompt');
        ylabel(handles(1, 1), 'Q prompt');

        %----- Nav bits ---------------------------------------------------
        plot  (handles(1, 2),trackResults.I_P');
        grid  (handles(1, 2));
        title (handles(1, 2), 'Bits of the navigation message');
        xlabel(handles(1, 2), 'Time (s)');
        axis  (handles(1, 2), 'tight');
        legend(handles(1, 2), 'I_{P}');

        %----- PLL discriminator unfiltered--------------------------------
        plot  (handles(2, 1),  ...
                              trackResults.pllDiscr, 'r');      

        grid  (handles(2, 1));
        axis  (handles(2, 1), 'tight');
        xlabel(handles(2, 1), 'Time (s)');
        ylabel(handles(2, 1), 'Amplitude');
        title (handles(2, 1), 'Raw PLL discriminator');


        %----- PLL discriminator filtered----------------------------------
        plot  (handles(3, 1),  ...
                              trackResults.pllDiscrFilt, 'b');      

        grid  (handles(3, 1));
        axis  (handles(3, 1), 'tight');
        xlabel(handles(3, 1), 'Time (s)');
        ylabel(handles(3, 1), 'Amplitude');
        title (handles(3, 1), 'Filtered PLL discriminator');
end
