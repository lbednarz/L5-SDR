function acqResults = acquisition(longSignal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end 
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the 
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
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
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Initialization =========================================================

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

% Create two 1msec vectors of data to correlate with and one with zero DC
signal1 = longSignal(1 : samplesPerCode);
signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);
% signal1 = longSignal(s*samplesPerCode : (s+1)*samplesPerCode-1);
% signal2 = longSignal((s+1)*samplesPerCode+1 : (s+2)*samplesPerCode);

% signal11 = longSignal(1 : 5*samplesPerCode);
% signal22 = longSignal(5*samplesPerCode+1 : 10*samplesPerCode);


signal0DC = longSignal - mean(longSignal);   %%Problems here....

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;
% phasePoints1 = (0 : (5*samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round(settings.acqSearchBand * 2) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable(settings);


%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);
% results1     = zeros(numberOfFrqBins, samplesPerCode);
CN1     = zeros(numberOfFrqBins, samplesPerCode);
% CN     = zeros(numberOfFrqBins, 2*samplesPerCode-1);


% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);

fprintf('(');

% Perform search for all listed PRN numbers ...
for PRN = settings.acqSatelliteList

%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));

    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins

        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               0.5e3 * (frqBinIndex - 1);

        
%         frqBins(frqBinIndex) = settings.IF +  20*(frqBinIndex + 1);
%                                (settings.acqSearchBand/2) * 1000 + ...
%                                0.5e3 * (frqBinIndex - 1);
% 
        %--- Generate local sine and cosine -------------------------------
        sigCarr = exp(1i*frqBins(frqBinIndex) * phasePoints);
%         sigCarr1 = exp(i*frqBins(frqBinIndex) * phasePoints1);
        
        %--- "Remove carrier" from the signal -----------------------------
        I1      = real(sigCarr .* signal1);
        Q1      = imag(sigCarr .* signal1);
        I2      = real(sigCarr .* signal2);
        Q2      = imag(sigCarr .* signal2);
        
%         I11      = real(sigCarr1 .* signal11);
%         Q11      = imag(sigCarr1 .* signal11);
%         I22      = real(sigCarr1 .* signal22);
%         Q22      = imag(sigCarr1 .* signal22);

        %--- Convert the baseband signal to frequency domain --------------
        IQfreqDom1 = fft(I1 + 1i*Q1);
        IQfreqDom2 = fft(I2 + 1i*Q2);
        
%         IQfreqDom11 = fft(I11 + j*Q11);
%         IQfreqDom22 = fft(I22 + j*Q22);

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        convCodeIQ1 = IQfreqDom1 .* caCodeFreqDom;
        convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;
        
%         convCodeIQ11 = IQfreqDom11 .* caCodeFreqDom;
%         convCodeIQ22 = IQfreqDom22 .* caCodeFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        r1 = (ifft(convCodeIQ1));
        r2 = (ifft(convCodeIQ2));
%         i1 = imag(ifft(convCodeIQ1));
%         i2 = imag(ifft(convCodeIQ2));
%         p1 = angle(ifft(convCodeIQ1));
%         p2 = angle(ifft(convCodeIQ2));
%         p1 = unwrap(p1,[],2);
%         p2 = unwrap(p2,[],2);
        acqRes1 = abs(ifft(convCodeIQ1)) .^ 2;
        acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
%         correlation_I1 = xcorr(I1,caCodesTable(PRN, :));
%         correlation_Q1 = xcorr(Q1,caCodesTable(PRN, :));
%         correlation1 = correlation_I1 + j*correlation_Q1;
% 
%         correlation_I2 = xcorr(I2,caCodesTable(PRN, :));
%         correlation_Q2 = xcorr(Q2,caCodesTable(PRN, :));
%         correlation2 = correlation_I2 + j*correlation_Q2;
        
%         acqRes11 = abs(ifft(convCodeIQ11)) .^ 2;
%         acqRes22 = abs(ifft(convCodeIQ22)) .^ 2;
%         acqRes1 = (r1+i1).^2;
%         acqRes2 = (r2+i2).^2;
%         acqRes1 = r1;
%         acqRes2 = r2;
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
%         if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
%             results1(frqBinIndex, :) = acqRes11;
            CN1(frqBinIndex, :) = r1;
%             CN(frqBinIndex, :) = correlation1;
%         else
%             results(frqBinIndex, :) = acqRes2;
% %             results1(frqBinIndex, :) = acqRes22;
%              CN1(frqBinIndex, :) = r2;
% %              CN(frqBinIndex, :) = correlation2;
%         end
        
%         results = abs(results);
    end % frqBinIndex = 1:numberOfFrqBins
% C = abs(CN); 
%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--- Find the correlation peak and the carrier frequency --------------
    [peakSize, frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [peakSize codePhase] = max(max(results));
    P_max = peakSize;
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = codePhase - samplesPerCodeChip;
    excludeRangeIndex2 = codePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    secondPeakSize = max(results(frequencyBinIndex, codePhaseRange));
    P_noise = mean(results(frequencyBinIndex, codePhaseRange));

    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(PRN) = peakSize/secondPeakSize;
    CP = results(frequencyBinIndex,:);
    % If the result is above threshold, then there is a signal ...
    if (peakSize/secondPeakSize) > settings.acqThreshold

%% Fine resolution frequency search =======================================
        
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('%02d ', PRN);
        
        %--- Generate 10msec long C/A codes sequence for given PRN --------
        caCode = generateCAcode(PRN);
        
        codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
                               (1/settings.codeFreqBasis));
                           
        longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
        CP = results(frequencyBinIndex,:);
        %--- Remove C/A code modulation from the original signal ----------
        % (Using detected C/A code phase)
        xCarrier = ...
            signal0DC(codePhase:(codePhase + 10*samplesPerCode-1)) ...
            .* longCaCode;
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency
        
        %--- Find the next highest power of two and increase by 8x --------
        fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
        
        %--- Compute the magnitude of the FFT, find maximum and the
        %associated carrier frequency 
        fftxc = abs(fft(xCarrier, fftNumPts)); 
        
        
        uniqFftPts = ceil((fftNumPts + 1) / 2);
        [fftMax, fftMaxIndex] = max(fftxc);
        fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
        if (fftMaxIndex > uniqFftPts) %and should validate using complex data
            if (rem(fftNumPts,2)==0)  %even number of points, so DC and Fs/2 computed
                fftFreqBinsRev=-fftFreqBins((uniqFftPts-1):-1:2);
                [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
                acqResults.carrFreq(PRN)  = -fftFreqBinsRev(fftMaxIndex);
            else  %odd points so only DC is not included
                fftFreqBinsRev=-fftFreqBins((uniqFftPts):-1:2);
                [fftMax, fftMaxIndex] = max(fftxc((uniqFftPts+1):length(fftxc)));
                acqResults.carrFreq(PRN)  = fftFreqBinsRev(fftMaxIndex);
            end
        else
            acqResults.carrFreq(PRN)  = (-1)^(settings.fileType-1)*fftFreqBins(fftMaxIndex);
        end
        
        acqResults.codePhase(PRN) = codePhase;
    
    else
        %--- No signal with this PRN --------------------------------------
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(')\n');
