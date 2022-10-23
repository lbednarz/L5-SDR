 function [trackResults, channel]= trackingL5asE1(fid, channel, settings)
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
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
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
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================

% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, settings.msToProcess);

% Freq of the C/A code:
trackResults.codeFreq       = inf(1, settings.msToProcess);

% Generate the Table of L5 code
[L5ICodesTable,L5QCodesTable] = makeL5Table(settings);


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

% Initialize filter parameters
for channelNr = 1:settings.numberOfChannels
    
        %FLL code tracking loop parameters
        Tracking.trackResults(channelNr).trackingweOldNco_1=0;
        Tracking.trackResults(channelNr).trackingweOldNco_2=0;
        Tracking.trackResults(channelNr).trackingweOldError_1=0;
        Tracking.trackResults(channelNr).trackingweOldError_2=0;

        Tracking.trackResults(channelNr).I_P_1=0;
        Tracking.trackResults(channelNr).Q_P_1=0;
        
        % PLL carrier/Costas loop parameters
        Tracking.trackResults(channelNr).trackingcarrOldNco_1=0;
        Tracking.trackResults(channelNr).trackingcarrOldNco_2=0;
        Tracking.trackResults(channelNr).trackingcarrOldError_1=0;
        Tracking.trackResults(channelNr).trackingcarrOldError_2=0;
        
        %DLL code tracking loop parameters
        Tracking.trackResults(channelNr).trackingcodeOldNco_1=0;
        Tracking.trackResults(channelNr).trackingcodeOldError_1=0;
   
end

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

trackTimes = settings.msToProcess;     % For GPS one C/A code is one ms

%--- DLL variables --------------------------------------------------------
% Define early-late offset (in chips)
earlyLateSpc = settings.dllCorrelatorSpacing;

% Summation interval (to say that the caculation is done every 1ms)
PDIcode = 0.001;

% Calculate filter coefficient values (calculate the loop filter coeficients)
[tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
                                    settings.dllDampingRatio, ...
                                    1.0);

%--- PLL variables --------------------------------------------------------
% Summation interval
PDIcarr = 0.001;

% Calculate filter coefficient values (calculate the loop filter coeficients)
[tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                    settings.pllDampingRatio, ...
                                    0.25);
hwb = waitbar(0,'Tracking...');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

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
            dataAdaptCoeff*(settings.skipNumberOfBytes + 2*(channel(channelNr).codePhase-1)), ...
            'bof');


        % Get a vector with the C/A code sampled 1x/chip
        L5dataCode  = L5ICodesTable(channel(channelNr).PRN, :);
        L5pilotCode = L5QCodesTable(channel(channelNr).PRN, :);
        % Then make it possible to do early and late versions
        L5dataCode  = [L5dataCode(end)  L5dataCode  L5dataCode(1)];
        L5pilotCode = [L5pilotCode(end) L5pilotCode L5pilotCode(1)];

        %--- Perform various initializations ------------------------------

        %initialize carrier setting
        Tracking.trackResults(channelNr).carrFreq = settings.IF;
        Tracking.trackResults(channelNr).carrierDoppler = channel(channelNr).acquiredFreq;
        Tracking.trackResults(channelNr).remCarrPhase = 0;
        
        %initialize code Freq setting
        Tracking.trackResults(channelNr).codeDoppler = Tracking.trackResults(channelNr).carrierDoppler * ...
             settings.codeFreqBasis / settings.carrFreqBasis; 
        Tracking.trackResults(channelNr).remCodePhase  = 0;
        Tracking.trackResults(channelNr).I_P_1=0;
        Tracking.trackResults(channelNr).Q_P_1=0;

        %=== Process the number of specified code periods =================
         for loopCnt =  1:trackTimes
            
            %FLL carrier tracking loop parameters
            weOldNco_1    = Tracking.trackResults(channelNr).trackingweOldNco_1;
            weOldNco_2    = Tracking.trackResults(channelNr).trackingweOldNco_2;
            weOldError_1  = Tracking.trackResults(channelNr).trackingweOldError_1;
            weOldError_2  = Tracking.trackResults(channelNr).trackingweOldError_2;

            % PLL carrier/Costas loop parameters
            carrOldNco_1   = Tracking.trackResults(channelNr).trackingcarrOldNco_1;
            carrOldNco_2   = Tracking.trackResults(channelNr).trackingcarrOldNco_2;
            carrOldError_1 = Tracking.trackResults(channelNr).trackingcarrOldError_1;
            carrOldError_2 = Tracking.trackResults(channelNr).trackingcarrOldError_2;

            %DLL code tracking loop parameters
            codeOldNco_1   = Tracking.trackResults(channelNr).trackingcodeOldNco_1;
            codeOldError_1 = Tracking.trackResults(channelNr).trackingcodeOldError_1;

            I_P_1=Tracking.trackResults(channelNr).I_P_1;
            Q_P_1=Tracking.trackResults(channelNr).Q_P_1;
             
%% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task. (IT IS THE BAR PROSSECING)
            if (rem(loopCnt, 50) == 0)
                try
                    waitbar(loopCnt/trackTimes, ...
                            hwb, ...
                            ['Tracking: Ch ', int2str(channelNr), ...
                            ' of ', int2str(settings.numberOfChannels), ...
                            '; PRN#', int2str(channel(channelNr).PRN), ...
                            '; Completed ',int2str(loopCnt), ...
                            ' of ', int2str(trackTimes), ' msec']);                       
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
            codeFrqNco = Tracking.trackResults(channelNr).codeDoppler + settings.codeFreqBasis;  %multiply by 12
            remCodePhase = Tracking.trackResults(channelNr).remCodePhase;
            codePhaseStep = codeFrqNco / settings.samplingFreq;
            
            % Number of sample in the code
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);
            
            % Read in the appropriate number of samples to process this
            % interation 
            [rawSignal, samplesRead] = fread(fid, ...
                dataAdaptCoeff*blksize, settings.dataType);
 
            rawSignal = rawSignal'; %Transpose the vector

            if (dataAdaptCoeff==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;  %Transpose vector
            end
            
            
            % If did not read in enough samples, then could be out of 
            % data - better exit 
            if (samplesRead ~= dataAdaptCoeff*blksize)
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
            earlyCode   = L5pilotCode(tcode2);
            
            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = L5pilotCode(tcode2);
            
            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                          codePhaseStep : ...
                          ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = L5pilotCode(tcode2);
            
            %Define index into prompt code vector for data channel
            promptCodeData =L5dataCode(tcode2);
            remCodePhase_New = (tcode(blksize) + codePhaseStep) - 10230.0;
            

%% Generate the carrier frequency to mix the signal to baseband -----------
            time    = (0:blksize) ./ settings.samplingFreq;
            carrFreq = settings.IF + Tracking.trackResults(channelNr).carrierDoppler;
            remCarrPhase = Tracking.trackResults(channelNr).remCarrPhase;
            
            % Get the argument to sin/cos functions
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase_New = rem(trigarg(blksize+1), (2 * pi));
            
            % Finally compute the signal to mix the collected data to
            % bandband
            carrsig = exp(1i .* trigarg(1:blksize));

%% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            
            data_I_P = sum(promptCodeData .* iBasebandSignal);
            data_Q_P = sum(promptCodeData .* qBasebandSignal);

            
%% Find FLL error and update carrier NCO ----------------------------------
                % FLL is not stable in this SDR....
                Pdot    = I_P_1*I_P+Q_P_1*Q_P;
                Pcross  = I_P_1*Q_P-Q_P_1*I_P;
    %             weError = Pcross*sign(Pdot)/(sqrt(I_P_1^2+Q_P_1^2)*sqrt(I_P^2+Q_P^2))/PDIcarr;
                weError = atan2(Pcross,Pdot)/PDIcarr;

                if strcmp(settings.TrackingSettingFLLUsed,'true')
                    weNco=FLLLoopFilter(PDIcarr,...
                        settings.FLLa2,settings.FLLBL2,...
                        weError,weOldError_1,weOldError_2,weOldNco_1,weOldNco_2);
                end
                
                
%% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P); %/(2*pi)
            
            % Implement carrier loop filter and generate NCO command
            carrNco=loopFilter2(PDIcarr,...
                    settings.PLLa3,settings.PLLb3,settings.PLLBL3,...
                    carrError,carrOldError_1,carrOldError_2,carrOldNco_1,carrOldNco_2);
                
            %% PLL+FLL carrier frequency controller
                if loopCnt <0
                     a=0;
                     b=1;
                else
                    a=1;
                    b=0;
                end         
                carrFrqNco=a*0.2*carrNco+b*0.07*weNco;
                carrierDoppler_new = channel(channelNr).acquiredFreq + carrFrqNco;

%% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            
            % Implement code loop filter and generate NCO command
            codeNco = codeOldNco_1 + (tau2code/tau1code) * ...
                (codeError - codeOldError_1) + codeError * (PDIcode/tau1code);

            
            % Modify code freq based on NCO command
            codeFrqNco = settings.codeFreqBasis - codeNco;
            codeDoppler_new = carrierDoppler_new  * 12 * settings.codeFreqBasis/ settings.carrFreqBasis + codeFrqNco ;

%% Record various measures to show in postprocessing ----------------------

            Tracking.trackResults(channelNr).I_P_1=I_P;
            Tracking.trackResults(channelNr).Q_P_1=Q_P;
            Tracking.trackResults(channelNr).codePhase=1;

            Tracking.trackResults(channelNr).carrFreq     = carrFreq;
            Tracking.trackResults(channelNr).remCarrPhase = remCarrPhase_New;
            Tracking.trackResults(channelNr).remCodePhase = remCodePhase_New;

            Tracking.trackResults(channelNr).carrierDoppler    = carrierDoppler_new;
            Tracking.trackResults(channelNr).codeDoppler       = codeDoppler_new;

            %FLL tracking loop parameters
            Tracking.trackResults(channelNr).trackingweOldNco_2 = Tracking.trackResults(channelNr).trackingweOldNco_1;
            Tracking.trackResults(channelNr).trackingweOldNco_1 = weNco;
            Tracking.trackResults(channelNr).trackingweOldError_2 = Tracking.trackResults(channelNr).trackingweOldError_1;
            Tracking.trackResults(channelNr).trackingweOldError_1 = weError;

            %DLL code tracking loop parameters
            Tracking.trackResults(channelNr).trackingcodeOldNco_1   = codeNco;
            Tracking.trackResults(channelNr).trackingcodeOldError_1 = codeError;

            %PLL carrier/Costas loop parameters
            Tracking.trackResults(channelNr).trackingcarrOldNco_2   = Tracking.trackResults(channelNr).trackingcarrOldNco_1;
            Tracking.trackResults(channelNr).trackingcarrOldNco_1   = carrNco;
            Tracking.trackResults(channelNr).trackingcarrOldError_2 = Tracking.trackResults(channelNr).trackingcarrOldError_1;
            Tracking.trackResults(channelNr).trackingcarrOldError_1 = carrError;

            % Record sample number (based on 8bit samples)
            trackResults(channelNr).absoluteSample(loopCnt) =(ftell(fid))/(2*dataAdaptCoeff);

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
            
            trackResults(channelNr).data_I_P(loopCnt) = data_I_P;
            trackResults(channelNr).data_Q_P(loopCnt) = data_Q_P;
            
        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;        
        
    end % if a PRN is assigned
end % for channelNr 

% Close the waitbar
close(hwb)
