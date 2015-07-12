clear all
KbName('UnifyKeyNames')
Screen('Preference', 'SkipSyncTests', 1);
input('hit enter to begin...  ');
%%%%%%%

output.Subject = 'test';

%%%%%%%
[keyboardIndices, productNames, ~] = GetKeyboardIndices;
% deviceString= 'Teensy Keyboard/Mouse';
deviceString = 'Apple Internal Keyboard / Trackpad'
for i=1:length(productNames)%for each possible device
    if strcmp(productNames{i},deviceString)%compare the name to the name you want
        deviceNumber=keyboardIndices(i);%grab the correct id, and exit loop
        break;
    end
end
if deviceNumber==0%%error checking
    error('No device by that name was detected');
end
triggerKey = 46;                % KbName('=+');
keylist = ones(1,256);          % keys for KbQueueCreate
keylist(triggerKey) = 0;        % dont want to record trigger
keyPressNumbers = {'82','81'};   % 82: up, 81: down

%%% Screen Parameters
w.whichScreen = 0;
w.ScreenWidth = 29;         % horizontal display size - 44.5;
w.ViewDistance = 57;        % in cm, ideal distance: 1 cm equals 1 visual degree (at 57 cm) - 98 at scanner
w.frameRate = 60;
w.ScreenSizePixels = Screen('Rect', w.whichScreen); %Scanner display = [0 0 1024 768];
w.VisAngle = (2*atan2(w.ScreenWidth/2, w.ViewDistance))*(180/pi); % Visual angle of the whole screen
stim.ppd = round(w.ScreenSizePixels(3)/w.VisAngle); % pixels per degree visual angle
sizepixel = (2*atan2((w.ScreenWidth/w.ScreenSizePixels(3))/2, w.ViewDistance))*(180/pi); % in visual degree
fNyquist = 0.5/sizepixel;

%%% STIMULUS PARAMETERS
stim.size = 15;         % in visual degree
stim.annulus = 3;
stim.fixation = 0.7;
stim.diameter = round(stim.size*stim.ppd); % in pixels
stim.annulus_diam = round(stim.annulus*stim.ppd);
stim.fix_diam = round(stim.fixation*stim.ppd);
stim.outer_fixation = round(1.4 * stim.fix_diam);

stim.rotAngle = 135;    %linspace(0,180,t.npresentation);
stim.filterwidth = 5;
stim.orient_filterWidthLow = stim.rotAngle - stim.filterwidth;
stim.orient_filterWidthHigh = stim.rotAngle + stim.filterwidth;
stim.Freq = [2 3];    % cpd
stim.ContrastNoise = .3; % max is 50% -> adding two signal together
stim.Grey = 128;
stim.AmplNoise = stim.ContrastNoise * (stim.Grey-1);

stim.phaseRad = 0;
stim.Ftarg = 4;  % cycles per degree
xysize = round(stim.diameter);
stim.xtrg =  0;   % relative to the center (in pixels) (positive:right)
stim.ytrgUp = -300;     % in pixels. positive (upward relative to the center)
stim.ytrgDown = 300;
stim.pUp = .5;
stim.target_alpha = pi/2;  % in radian
stim.MaxContrTarg = .9 - stim.ContrastNoise;

f = freqspace(stim.diameter);
Freq = stim.Freq .* stim.size;

CenterX = w.ScreenSizePixels(3)/2; CenterY = w.ScreenSizePixels(4)/2;

%%% Timing parameters:
t.TheDate = datestr(now,'yymmdd');  %Collect todays date
t.TimeStamp = datestr(now,'HHMM');  %Timestamp for saving out a uniquely named datafile (so you will never accidentally overwrite stuff)
t.stimpresentation = 2;             % stimulus presentation time
t.responseDur = 1;                  % response time
t.npresentation = 10;                % number of stimulus presentation in each block
t.nframesStim = t.stimpresentation*w.frameRate;
t.nframesResp = t.responseDur*w.frameRate;
t.init = 1;
StartTimes = t.init:(t.stimpresentation+t.responseDur):t.npresentation*(t.stimpresentation+t.responseDur);
EndTimes = StartTimes + t.stimpresentation;

%%% CREATE ORIENTATION FILTERED NOISE (vertical)
% Setup Filter (spatial freq and orientation)
[x,y] = meshgrid((-xysize/2):(xysize/2)-1, (-xysize/2):(xysize/2)-1);
freqFilter = Bandpass2(xysize, f(Freq(1)), f(Freq(2)));
orientFilter = OrientationBandpass(xysize,stim.orient_filterWidthLow, stim.orient_filterWidthHigh);
Filter = freqFilter .* orientFilter; Filter(Filter> 0) = 1;
h = fspecial('disk', 3); Filter = conv2(Filter, h, 'same');
Filter = Filter./max(Filter(:));
% makes sure dc component is not filtered out:
temp_filter = fftshift(Filter); temp_filter(1,1) = 1; Filter = ifftshift(temp_filter);

% Gaussian Mask
gaussian_std = round(stim.ppd*3);
eccen = sqrt((x).^2+(y).^2); 	% calculate eccentricity of each point in grid relative to center of 2D image
Gaussian = zeros(xysize); Gaussian(eccen <= (xysize/2)-stim.ppd/2) = 1;
Gaussian = conv2(Gaussian, fspecial('gaussian', stim.ppd, stim.ppd), 'same');

% Gaussian Annulus
[X,Y] = meshgrid(1:xysize,1:xysize);
Annulus = zeros(xysize);
r_eccen = sqrt((X-xysize/2).^2+(Y-xysize/2).^2); 	% calculate eccentricity of each point in grid relative to center
Annulus(r_eccen > stim.annulus_diam/2) = 1;
Annulus = conv2(Annulus, fspecial('gaussian', 30, 10), 'same');
Annulus(r_eccen > stim.annulus_diam) = 1;


for n = 1:t.npresentation
    noise_tmp = 2*((rand(xysize, xysize)-.5 ));
    fn = fftshift(fft2(noise_tmp));
    filteredNoise = real(ifft2(ifftshift(Filter.*fn)));
    noise = filteredNoise./max(max(abs(filteredNoise))); % scales from [-1 1], contrast can be adjusted later
    noiseField{n} = noise .* Gaussian .* Annulus;
end

% Gaussian Target
Pti = 1:xysize;
[Xt,Yt] = meshgrid(Pti,Pti);
SinTarg = sin( 2.*pi.*stim.Ftarg.*sizepixel.*cos(stim.target_alpha).*Xt + 2.*pi.*stim.Ftarg.*sizepixel.*sin(stim.target_alpha).*Yt + stim.phaseRad);

gausTarg_std = round(stim.ppd/2);
eccen_targ = (x-stim.xtrg).^2+(y-stim.ytrgUp).^2;
GaussianTargUp = exp(-eccen_targ./(2.*gausTarg_std^2)); 	% calculate eccentricity of each point in grid relative to center of 2D image

clear eccen_targ
eccen_targ = (x-stim.xtrg).^2+(y-stim.ytrgDown).^2;
GaussianTargDown = exp(-eccen_targ./(2.*gausTarg_std^2)); 	% calculate eccentricity of each point in grid relative to center of 2D image

TargetUp  = GaussianTargUp.*SinTarg;
TargetDown  = GaussianTargDown.*SinTarg;

%%% Create Trial Events
%
for ii=1:t.npresentation
    if rand<stim.pUp
        stimUpTrials(ii) = 1;
    end
end
%}
%stimUpTrials = ones(1,t.npresentation);
%%%%%% Start the experiment %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize the run:

%
noiseField = Shuffle(cellfun(@flipud, noiseField, 'UniformOutput', false));
% initializing the QUEST procedure:
quest.pThreshold=0.95;
quest.beta=1;quest.delta=0.01;quest.gamma=0.5;quest.tGuess = -1; quest.tGuessSd = 2;
grain = [];range = log10(stim.MaxContrTarg * stim.Grey); plotIt = 1;
q=QuestCreate(quest.tGuess,quest.tGuessSd,quest.pThreshold,quest.beta,quest.delta,quest.gamma,grain, range,plotIt);
response = [1 1 1 1 1 1 1 1];
%%% WINDOW SETUP
AssertOpenGL;
[window rect] = Screen('OpenWindow',w.whichScreen, stim.Grey); %,[0 0 400 400]);
HideCursor;
%%% MAKE COLOR LOOKUP TABLE AND APPLY GAMMA CORRECTION
OriginalCLUT = Screen('LoadCLUT', window);
yellow = [255 255 0]; red = [255 0 0]; green = [0 255 0]; blue = [0 0 255];
white = WhiteIndex(w.whichScreen);
black = BlackIndex(w.whichScreen);
%%%
Screen('TextStyle', window, 1);
Screen('TextSize', window, 16);
bbox = Screen('TextBounds', window, 'X');
newRect = CenterRectOnPoint(bbox, CenterX, CenterY);
tx = newRect(1); ty = newRect(2);
Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));
Screen('DrawText', window, 'Start', tx-10, ty + 40, 255, 0);
Screen('Flip', window);

GetClicks;

PsychHID('KbQueueCreate', deviceNumber, keylist);

StartTimeExp = GetSecs;
threshUpdate = nan(1,t.npresentation);
%%%
for n=1:t.npresentation
    tempStim = noiseField{n};
    tempNoise = tempStim./max(tempStim(:)).* stim.AmplNoise ;
    
    %%%  setting the next QUEST test threshold (choose one of these three methods):
    tTest=QuestQuantile(q);	% Recommended by Pelli (1987), and still our favorite.
    % 	tTest=QuestMean(q);		% Recommended by King-Smith et al. (1994)
    % 	tTest=QuestMode(q);		% Recommended by Watson & Pelli (1983)
    
    %%%
    stim.AmplTarg(n) = 10^(q.xThreshold);
    %%%
    if ~stimUpTrials(n)
        TargetDown = TargetDown./max(abs(TargetDown(:))).*stim.AmplTarg(n);
        noisyStim = tempNoise + TargetDown + stim.Grey;
        UpOn = 0;
    else
        TargetUp = TargetUp./max(abs(TargetUp(:))).*stim.AmplTarg(n);
        noisyStim = tempNoise + TargetUp + stim.Grey;
        UpOn = 1;
    end
    Stimulus = Screen('MakeTexture', window, noisyStim);
    for kk=1:round(t.nframesStim)
        if GetSecs-StartTimeExp> EndTimes(n)
            break
        else
            Screen('DrawTexture', window,Stimulus);
            Screen('Flip', window, StartTimeExp+StartTimes(n), [], [], 0);
        end
    end
    KbQueueFlush();
    Screen('Close',Stimulus)
    %%% finding the response:
    PsychHID('KbQueueStart', deviceNumber);
    responsewindowTrig_RSVP = 1;
    for kk=1:round(t.nframesResp)
        GetSecs-StartTimeExp
        t.responseDur+EndTimes(n)
        if GetSecs-StartTimeExp>t.responseDur+EndTimes(n)
            break
        else
            Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
            Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));
            Screen('DrawText', window, 'Up or Down?', tx-10, ty + 40, 255, 0);
            Screen('Flip', window);
        end
    end
    [pressed, firstpress] = PsychHID('KbQueueCheck', deviceNumber);
    whichkeys = find(firstpress);
    if ((strcmp(num2str(whichkeys), keyPressNumbers{1}) && UpOn) || (strcmp(num2str(whichkeys), keyPressNumbers{2}) && ~UpOn))
        RSVP_rightwrong(n) = 1;
    else
        RSVP_rightwrong(n) = 0;
    end
    %%% update QUEST procedure:
    q=QuestUpdate(q,tTest,RSVP_rightwrong(n));
end

EndTime = GetSecs;
quest.EstThresh = QuestMean(q);		% or QuestMode or QuestQuantile
quest.EstSd = QuestSd(q);
fprintf('Final threshold estimate (mean+-sd) is %.2f +- %.2f\n',quest.EstThresh,quest.EstSd);


%{
%%% make & output
output.RSVP_rightwrong = RSVP_rightwrong;
output.mean_RSVP_rightwrong = mean(RSVP_rightwrong);
TheData{runnumber}.stim = stim;
TheData{runnumber}.t = t;
TheData{runnumber}.w = w;
TheData{runnumber}.mylog.stimtimes_s = {output.trueOn_Time};
TheData{runnumber}.mylog.durationss_s = {output.trueOn_duration};
TheData{runnumber}.mylog.designMatrix = output.designMatrix;
TheData{runnumber}.output = output;
eval(['save ', 'Data_tunedSuppression_', output.Subject, '.mat TheData']);
Screen('FillOval', window, white, CenterRectOnPoint([0 0 stim.outer_fixation stim.outer_fixation], CenterX, CenterY));
Screen('FillOval', window, stim.Grey, CenterRectOnPoint([0 0 stim.fix_diam stim.fix_diam], CenterX, CenterY));
Screen('DrawText', window, 'DONE.', CenterX-20, CenterY-stim.annulus_diam/2);
Screen('Flip', window);

WaitSecs(3)
%}
ShowCursor;
Screen('CloseAll')