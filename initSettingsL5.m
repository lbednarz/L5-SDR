function settings = initSettingsL5()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".  
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 

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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided

settings.msToProcess        = 36000; %420*1000;        %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels   = 8;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only.
settings.samplesToskip   = 0;
settings.msToskip        = 15;
% 400267;
% settings.prn             = 10;
% settings.msToskip        = 811;
settings.skipNumberOfBytes     = settings.msToskip*16*25000/8;

% settings.skipNumberOfBytes     = 76029*20000*16/8;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode
settings.fileName           = 'C:\Users\Sahil\Desktop\Sterenn\Matlab\test_L5_20220706_RX2.bin';
% 'C:\Users\Sahil\Desktop\Sterenn\Matlab\B200Ex_1_Fif20k_Fs18M.bin';
% 'C:\Users\Sahil\Desktop\Sterenn\Matlab\Test_IQ_Fc1157_Fs78MHz_mat.bin';
% 'C:\Users\Sahil\Desktop\Sterenn\Matlab\Test_IQ_Fc1153_Fs78MHz.bin';
% 'C:\Users\Sahil\Desktop\Sterenn\Matlab\Test_IQ_Fc1157_Fs78MHz.bin';
% 'C:\Users\Sahil\Desktop\Sterenn\Matlab\test_L5_20220706_RX2.bin';
% 'C:\Users\Sahil\Desktop\Sterenn\Matlab\test_L5_20220722_Rx2.bin';
% 'C:\Users\Sahil\Desktop\GNSS\cleanStatic_TEX.bin'; 

   
% Data type used to store one sample
settings.dataType           = 'int16';


% File Types
%1 - 16 bit real samples S0,S1,S2,...
%2 - 16 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...                      
settings.fileType           = 2;

% Intermediate, sampling and code frequencies
settings.IF                 = 0e3;      %[Hz]
settings.samplingFreq       = 25e6;       %[Hz], respect the Nyquist-Shannon theorem (saw 20.48e6 Hz in comparaison L1/L5)
settings.codeFreqBasis      = 10.23e6;    %[Hz]
settings.carrFreqBasis      = 1176.45e6;  %[Hz]
settings.FreqSec            = 1000;       %[Hz]

% Define number of chips in a code period
settings.codeLength         = 10230;

%% Acquisition settings ===================================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
settings.acqSatelliteList   = 1:32;         %[PRN numbers]
% Band around IF to search for satellite signal. Depends on max Doppler
settings.acqSearchBand      = 14;           %[kHz] to change but I don't know what to put
% Threshold for the signal presence decision rule
settings.acqThreshold       = 2.5;

%% Tracking loops settings ================================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;       %[Hz]
settings.dllCorrelatorSpacing    = 0.5;     %[chips]
settings.dllVeryEarlyLateSpc     = 0.75;    %[chips]

% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 15;      %[Hz]
settings.PLLa3                   = 1.414;
settings.PLLb3                   = 2.4;
settings.PLLBL3                  = 15;

% FLL carrier tracking loop parameters
settings.TrackingSettingFLLUsed  =  'true';
settings.FLLa2                   =   1.414;
settings.FLLBL2                  =       8;      %[Hz]

%% Navigation solution settings ===========================================

% Period for calculating pseudoranges and position
settings.navSolPeriod       = 1;          %[ms]

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 10;           %[degrees 0 - 90] 
% Enable/dissable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off
                                            % 1 - On

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 1;            % 0 - Off
                                            % 1 - On
%% Nav bits estimation

settings.bitEst = 2; %
%% Constants ==============================================================

settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time
