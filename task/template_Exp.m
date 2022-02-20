% =======================
% function[] = BKH_Bubbles(run_mode)
% run_mode = 1: experiment
% run_mode = 2: two practice block with 10 trials each
% run_mode = 3: demo w/o responses
%
% GV 14.02.2021
% GV 17.02.2021
% GV 08.06.2021
% GV 01.07.2021
% =======================

function[] = BKH_Bubbles()
Screen('Preference', 'SkipSyncTests', 1);

%demo         = 1;    % runs demo without participant numbers etc
run_mode = 1;
if ~ismember(run_mode, [1 2 3])
    error('Wrong runmode argument\n'); end

% set up paths, responses, monitor, ...
[vp, condition_code, facedir, expdir, MonitorSelection, response_mapping, instruct] = get_experimentInfo(run_mode);
[TastenVector, ~, response_keys] = prepare_responseKeys(response_mapping); % manual responses
subdir    = check_dir([expdir, vp, '\']); % check for (and possibly create) subject data directory
lastblock = check_files(subdir, condition_code); % check for previous data files
MonitorSpecs = getMonitorSpecs(MonitorSelection); % subfunction, gets specification for monitor

% Bubble properties
num_cycles   = 2;  % number of cycles that the bubble patches will contain
Initnbubbles = 120; % start number of bubbles at scale with highest SF
sd_Gauss     = num_cycles / 2.355; % as a fraction of wavelength of cycles per face, cfp
                     % FWHM is 2.355*sd_Gauss
                     % 2 cpf at FWHM is sd = 2/2.355 = 0.849
                     % 3cpf at FWHM is 1.2739
                     % 2.5 cpf at FWHM is 1.0616
                     
% Quest procedure
tGuess     = log10(Initnbubbles); % start with this number of bubbles
tGuessSd   = 0.1;
pThreshold = 0.75;
beta = 3.5; delta = 0.01; gamma = 0.5;

% Display Options
faceSize     = 8;  % display size of face, in cycles per deg visual angle
use_lowpassresidual = 0; % 0 = none, 1 = individual, 2 = common


% read stimuli, construct gabor patches for scales
[stim] = get_stims([facedir, 'p5_struct_npic_470x349.mat'], condition_code); % see subfunction
%stim
[patches] = get_patches(stim.facedims, stim.mids, num_cycles, sd_Gauss);
stimSize = faceSize * (stim.picdims(1) / stim.facedims(1));
DestRectSize = round([1 1 stimSize * MonitorSpecs.PixelsPerDegree stimSize * MonitorSpecs.PixelsPerDegree]); %left, top, right, bottom


% for QUEST procedur; see QuestDemo.m (Variable naming is similar to QuestDemo.m)
if ~isempty(lastblock)
   load(lastblock);
   q = Bubbles.q; clear Bubbles;
else
    q = QuestCreate(tGuess, tGuessSd, pThreshold, beta, delta, gamma);
    q.normalizePdf = 1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
end


try
    Priority(1);
    HideCursor;
    [win, ~] = Screen('OpenWindow', MonitorSpecs.ScreenNumber, 127); 

    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [width, height] = Screen('WindowSize', win);
    hz = Screen('NominalFrameRate', win);
    Screen('Flip', win);
    Screen('TextSize', win, 25); %fontsize 34
    Screen('TextFont', win, 'Helvetica');

    dRect = CenterRectOnPoint(DestRectSize, width/2, height/2);    

    % Construct and start response queue
    KbQueueCreate([],TastenVector);
    KbQueueStart;

    
    blocknum = 0;
    instTexture = Screen('MakeTexture', win, instruct);
    escpress = show_instruction(win, instTexture, run_mode, blocknum);

    % block loop
    while ~escpress
       
       blocknum = blocknum + 1; 
       tstring = datestr(clock,30);  
       outfilename = [vp, '_', tstring, '.mat'];
       outm = []; Bubbles = []; 
       WaitSecs(1); 

       indx = Shuffle(1:numel(stim.names)); 
            if run_mode == 2
                indx = indx(1:10); end
     % trial loop
        for trial = 1:numel(indx)
        tbubbles  = QuestQuantile(q); % QUEST
        n_bubbles = round(10 ^ tbubbles);
                        
        [a_planes, b_centers, b_dims, f_coords] = get_alphaplanes(patches, stim.npic(indx(trial),:), n_bubbles, stim.facedims);
        [t_stim, ~] = sum_to_target(a_planes, stim.npic(indx(trial),:), use_lowpassresidual, stim.lp_avg);
        
        Texture = Screen('MakeTexture', win, t_stim);
        Screen('DrawTexture', win, Texture,[], dRect);
        FlipTime = Screen('Flip', win);
        
        %[stim.cond(indx(trial)), n_bubbles]
        %indx(trial)
        %stim.cond(indx(trial))

        [rkey_code, RT, IsResponseCorrect] = get_behavioralresponse(FlipTime, response_keys(stim.cond(indx(trial))), run_mode, n_bubbles, Initnbubbles);
        q = QuestUpdate(q, tbubbles, IsResponseCorrect);
        
        WaitSecs(0.2);
        Screen('Close', Texture);
        Screen('Flip', win);
        WaitSecs(1);
   
        %[stim.cond(indx(trial)), n_bubbles, IsResponseCorrect]
        outmat(trial,:) = [trial, indx(trial), rkey_code, RT, IsResponseCorrect, n_bubbles, ...
            stim.sex(indx(trial)), stim.emot(indx(trial)), ...
            stim.cond(indx(trial))];
        
        b_centers_cell{trial,:} = b_centers;
        b_dims_cell{trial,:}    = b_dims;
        f_coords_cell{trial,:}  = f_coords;
        
    end
        
    escpress = show_instruction(win, instTexture, run_mode, blocknum);
        
    Bubbles.outmat     = outmat;
    Bubbles.b_centers  = b_centers_cell;
    Bubbles.b_dims     = b_dims_cell;
    Bubbles.f_coords   = f_coords_cell;
    Bubbles.num_cycles = num_cycles;
    Bubbles.picdims    = stim.picdims;
    Bubbles.facedims   = stim.facedims;
    Bubbles.sd_Gauss   = sd_Gauss;
    Bubbles.vp         = vp;
    Bubbles.q          = q;
    Bubbles.date       = tstring;
    Bubbles.block      = blocknum;
    Bubbles.condition  = stim.stimstr;
    Bubbles.width      = width;
    Bubbles.height     = height;
    Bubbles.Hz         = hz;
    Bubbles.stmfile    = 'p5_struct_npic_470x349';
    
    save([expdir, vp, '\', outfilename], 'Bubbles');

    end % while
catch
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end
Screen('CloseAll');
end

%% =================== subfunctions =========

function [escpress] = show_instruction(windowPointer, instTexture, run_mode, blocknum)
       if run_mode == 2 && blocknum == 2
        escpress = 1;
       elseif run_mode == 3
       Screen('Flip', windowPointer);
       Screen('DrawTexture', windowPointer, instTexture);
       Screen('Flip', windowPointer); 
       WaitSecs(1); % stimulus presentation time is until response
       escpress = 0; 
       else 
       Screen('Flip', windowPointer);
       Screen('DrawTexture', windowPointer, instTexture);
       Screen('Flip', windowPointer); 
       KbQueueFlush();
       [~, keyCode, ~] = KbStrokeWait();
       escpress = (find(keyCode) == KbName('ESCAPE'));
       end
       WaitSecs(0.2);
       Screen('Flip', windowPointer);
end


function [rkey_code, RT, IsResponseCorrect] = get_behavioralresponse(FlipTime, ActuallyCorrect, run_mode, n_bubbles, Initnbubbles)


if run_mode == 3
   n_bubbles
   WaitSecs(1); % stimulus presentation time is until response
   rprob = 1./(1+exp(-0.1*(n_bubbles-Initnbubbles))); % log function with mean 120
   IsResponseCorrect = binornd(1, rprob);
   rkey_code = NaN;
   RT = NaN;
   [~, timevec] = KbQueueCheck;
   GetBehavioralPerformance(timevec); % for ESC press
else
    KbQueueFlush(); pressed = 0; % flush queue
    while ~pressed           
         [pressed, timevec] = KbQueueCheck; % no time limit for response, see Gosslyn & Schyns 2001 
    end
[KEYvec, RTvec] = GetBehavioralPerformance(timevec); % get RT; see subfunction
rkey_code = KEYvec(1);
RT = RTvec(1) - FlipTime;
IsResponseCorrect = rkey_code == ActuallyCorrect;
end
end

function [lastblock] = check_files(subdir, condition_code)
if condition_code == 1
   previousfiles = dir([subdir,'*NE-HA.mat']);
elseif condition_code == 2
   previousfiles = dir([subdir,'*NE-SA.mat']);
end

if ~ isempty(previousfiles)
tmpf = sort({previousfiles.name});
%howMany 
lastblock = tmpf{end};
else
 lastblock = [];   
end
end


function [subdir] = check_dir(subdir)
if ~ (exist(subdir, 'dir') == 7)
    mkdir(subdir); end
end


function [stim] = get_stims(stimFileName, condition_code) 
    tmpname = load(stimFileName);
    npic = tmpname.struct_npic;
    % neutral-happy
    pic_numbers = npic.picnum.h;
    NE_HA.npic  = npic.npic(pic_numbers,:);
    NE_HA.names = npic.names(pic_numbers);
    NE_HA.lp_avg  = npic.NE_HA_lpavg;
    tmp        = char(NE_HA.names);
    for j = 1:size(tmp, 1)
        tmp2 = deblank(tmp(j,:));
    NE_HA.sex(j)   = (tmp2(end-14) == 'F') + 1; % 1 = male, 2 = female
    NE_HA.emot(j)  = (tmp2(end-7) == 'S') + 1; % 1 = happy, 2 = sad
    NE_HA.cond(j)  = str2num(tmp2(end-4));     % 1 = neutral, 2 = emotion
    end
    %NE_HA.emot(NE_HA.cond == 1) = NaN;
    NE_HA.num   = pic_numbers'; clear tmp tmp2   
    % neutral-sad
    pic_numbers = npic.picnum.s;
    NE_SA.npic  = npic.npic(pic_numbers,:);
    NE_SA.names = npic.names(pic_numbers);
    NE_SA.lp_avg    = npic.NE_SA_lpavg;
    tmp        = char(NE_SA.names);
    for j = 1:size(tmp, 1)
        tmp2 = deblank(tmp(j,:));
    NE_SA.sex(j)   = (tmp2(end-14) == 'F') + 1; % 1 = male, 2 = female
    NE_SA.emot(j)  = (tmp2(end-7) == 'S') + 1; % 1 = happy, 2 = sad
    NE_SA.cond(j)  = str2num(tmp2(end-4));     % 1 = neutral, 2 = emotion
    end
    %NE_SA.emot(NE_SA.cond == 1) = NaN;
    NE_SA.num   = pic_numbers'; clear tmp2
    
    % return values
    if condition_code == 1
        stim    = NE_HA;
        stimstr = 'NE_HA';
    elseif condition_code == 2
        stim = NE_SA;
        stimstr = 'NE_SA';
    else
        error('\nWrong condition code!');
    end
    
%     picdims = npic.picdims;
%     facedims = npic.facedims;
%     mids    = npic.mids;
    
    stim.picdims  = npic.picdims;
    stim.facedims = npic.facedims;
    stim.mids     = npic.mids;
    stim.stimstr  = stimstr;
end



function [KEYvec, RTvec] = GetBehavioralPerformance(TimeVector)
indi=find(TimeVector); % wo ungleich 0, also responses
[~, i]=sort(TimeVector(indi)); % sortiere nach Zeit
KEYvec = indi(i); % key, aufsteigend nach Zeit
RTvec  = TimeVector(indi(i)); % rt, aufsteigend nach Zeit, relativ zum letzte KbEventFlush
% stop experiment if ESC is pressed
    if ismember(KbName('ESCAPE'), KEYvec) %ESC-taste
       sca; end
end

% function [gausswin2d] = NewGetGausswin(picsize, sdpix);
% %picsize=150;
% %sdpix   = 5.22% PixelsPerDegree * sdGauss;
% k = (sdpix/picsize) *2; 
% w1 = window(@gausswin,picsize,1/k); %previously k = 0.35
% gausswin2d = w1*w1';
% return
% end


function [MonitorSpecs] = getMonitorSpecs(MonitorSelection)
switch MonitorSelection
    case 1 % Dell-Monitor EEG-Labor
    MonitorSpecs.width     = 530; % in mm
    MonitorSpecs.height    = 295; % in mm
    MonitorSpecs.distance  = 600;%810; % in mm
    MonitorSpecs.xResolution = 2560; % x resolution
    MonitorSpecs.yResolution = 1440;  % y resolution
    MonitorSpecs.refresh     = 120;  % refresh rate
    MonitorSpecs.PixelsPerDegree = round((2*(tan(2*pi/360*(1/2))*MonitorSpecs.distance))/(MonitorSpecs.width / MonitorSpecs.xResolution)); % pixel per degree of visual angle
    case 2 % alter Belinea-Monitor Büro
    MonitorSpecs.width     = 390; % in mm
    MonitorSpecs.height    = 290; % in mm
    MonitorSpecs.distance  = 600; % in mm
    MonitorSpecs.xResolution = 1280; % x resolution
    MonitorSpecs.yResolution = 1024;  % y resolution
    MonitorSpecs.refresh     = 60;  % refresh rate
    MonitorSpecs.PixelsPerDegree = round((2*(tan(2*pi/360*(1/2))*MonitorSpecs.distance))/(MonitorSpecs.width / MonitorSpecs.xResolution)); % pixel per degree of visual angle
    case 3 % neuer Monitor Büro, Dell P2210
    MonitorSpecs.width     = 470; % in mm
    MonitorSpecs.height    = 295; % in mm
    MonitorSpecs.distance  = 600; % in mm
    MonitorSpecs.xResolution = 1680; % x resolution
    MonitorSpecs.yResolution = 1050;  % y resolution
    MonitorSpecs.refresh     = 60;  % refresh rate
    MonitorSpecs.PixelsPerDegree = round((2*(tan(2*pi/360*(1/2))*MonitorSpecs.distance))/(MonitorSpecs.width / MonitorSpecs.xResolution)); % pixel per degree of visual angle
    MonitorSpecs.ScreenNumber    = 0;
    case 4 % Monitor Home Office, LG 24BK750Y-B, 24-inch
    MonitorSpecs.width     = 530; % in mm
    MonitorSpecs.height    = 295; % in mm
    MonitorSpecs.distance  = 600; % in mm
    MonitorSpecs.xResolution = 1920; % x resolution
    MonitorSpecs.yResolution = 1080;  % y resolution
    MonitorSpecs.refresh     = 60;  % refresh rate
    MonitorSpecs.PixelsPerDegree = round((2*(tan(2*pi/360*(1/2))*MonitorSpecs.distance))/(MonitorSpecs.width / MonitorSpecs.xResolution)); % pixel per degree of visual angle
    MonitorSpecs.ScreenNumber    = 1;
    case 5 % Alex Notebook HP Elitebook 850 6G
    MonitorSpecs.width     = 345; % in mm
    MonitorSpecs.height    = 195; % in mm
    MonitorSpecs.distance  = 600; % in mm
    MonitorSpecs.xResolution = 1920; % x resolution
    MonitorSpecs.yResolution = 1080;  % y resolution
    MonitorSpecs.refresh     = 60;  % refresh rate
    MonitorSpecs.PixelsPerDegree = round((2*(tan(2*pi/360*(1/2))*MonitorSpecs.distance))/(MonitorSpecs.width / MonitorSpecs.xResolution)); % pixel per degree of visual angle
    MonitorSpecs.ScreenNumber    = 0;
end
end

function [vp, condition_code, facedir, expdir, MonSelection, response_mapping, instruct] = get_experimentInfo(run_mode)
if run_mode < 3
  vp = input('\nParticipant (three characters, e.g. S01)? ', 's');
    if length(vp)~=3 
       error ('Wrong input!'); end
  condition_code = str2num(input('\nCondition (1 = NE-HA, 2 = NE-SA)? ', 's'));  
    if ~ismember(condition_code, [1, 2])
        error('\nWrong condition code!'); end
  response_mapping = str2num(input('\nResponse mapping?\n1: y = neutral, m = emotion\n2: y = emotion, m = neutral\n', 's'));    
      if ~ismember(response_mapping, [1, 2])
        error('\nWrong condition code!'); end
else
      vp = 'tmp';
      condition_code = 1;
      response_mapping = 1;
end

if strcmp(getenv('COMPUTERNAME') , 'DESKTOP-MDU4B4V')
  facedir  = 'C:\Users\Gregor\Filr\Meine Dateien\Bubbles_BKH\';
  expdir   = 'C:\Users\Gregor\Filr\Meine Dateien\Bubbles_BKH\data\';
  MonSelection = 4;
elseif strcmp(getenv('COMPUTERNAME') , 'PC2022001273')
  facedir  = 'C:\Users\Alex\Meine Dateien\Bubbles_BKH\';
  expdir   = 'C:\Users\Alex\Meine Dateien\Bubbles_BKH\data\';
  MonSelection = 5;
end

switch response_mapping
    case  1
        instruct = imread([facedir, 'instruction_mapping1.png']);
    case  2
        instruct = imread([facedir, 'instruction_mapping2.png']);
end
instruct = 255-instruct; instruct(instruct == 255) = 127;
end

function[patch] = get_patches(picdims, mids, num_cycles, sd_Gauss)
 onecycle  = picdims(2) ./ fliplr(mids);
 allcycles = onecycle * num_cycles;

 sd_pic = onecycle .* sd_Gauss;
 picsize  = onecycle .* num_cycles * 2; % large pics so that gauss comes to zero
 k = (sd_pic./picsize) *2; 
        
 for md = 1:numel(mids)
    w1 = window(@gausswin, round(picsize(md)),1/k(md)); %previously k = 0.35
    patch{md} = w1*w1';
 end
end


function [t_stim, t_planes] = sum_to_target(a_planes, npic, use_lowpassresidual, lp_avg)
t_stim = zeros(size(npic{1}));
for m = 1:numel(a_planes)
t_planes{m} = npic{m}.*a_planes{m};
t_stim = t_stim + t_planes{m};
end
if use_lowpassresidual == 0
    t_stim = t_stim + 127;
elseif  use_lowpassresidual == 1
    t_stim = t_stim + npic{length(npic)};
elseif use_lowpassresidual == 2
    t_stim = t_stim + lp_avg;
%     tdims = size(common_lp);
%     dims  = size(t_stim);
%     margins = round(tdims/2 - dims/2);
%     common_lp(margins(1)+1 : (margins(1)+dims(1)), margins(2)+1 : (margins(2)+dims(2))) = ...
%         common_lp(margins(1)+1 : (margins(1)+dims(1)), margins(2)+1 : (margins(2)+dims(2))) + t_stim;
%     t_stim = common_lp;
end
end

function [a_planes, b_centers, b_dims, f_coords] = get_alphaplanes(patches, npics, nbubbles, t_area)
max_alpha = max(patches{1}(:));
for k = 1:(numel(npics)-1)
numbubbles = round((nbubbles / (4 ^ (k-1))));
[bubble_center, bubble_dims, face_coords] = prepare_alpha(patches{k}, npics{k}, numbubbles, t_area);
[alpha_plane] = add_alphaplane(patches{k}, npics{k}, bubble_center, bubble_dims, face_coords);
alpha_plane(alpha_plane > max_alpha) = max_alpha;
a_planes{k}  = alpha_plane;
b_centers{k} = bubble_center;
b_dims{k}    = bubble_dims;
f_coords{k}  = face_coords;
end
end
% 
% function [a_planes, b_centers, b_dims, f_coords] = get_alphaplanes(patches, npics, nbubbles, bubble_center, bubble_dims, face_coords);
% max_alpha = max(patches{1}(:));
% for k = 1:(numel(npics)-1)
% numbubbles = round((nbubbles / (4 ^ (k-1))));
% [bubble_center, bubble_dims, face_coords] = prepare_alpha(patches{k}, npics{k}, numbubbles);
% [alpha_plane] = add_alphaplane(patches{k}, npics{k}, bubble_center, bubble_dims, face_coords);
% alpha_plane(alpha_plane > max_alpha) = max_alpha;
% a_planes{k}  = alpha_plane;
% b_centers{k} = bubble_center;
% b_dims{k}    = bubble_dims;
% f_coords{k}  = face_coords;
% end
% end


function [bubble_center, bubble_dims, face_coords] = prepare_alpha(patch, face, n_bubbles, t_area)
% generate random locations for bubbles (center)
if any(t_area > size(face))
    error('impossible target indices'); end
offset = round((size(face) - t_area)/2);

if n_bubbles > 0
target_center = randperm(prod(t_area), n_bubbles);

    for bubble = 1:n_bubbles
    % location in face where bubble is to be centered
    [trow, tcol]  = ind2sub(t_area, target_center(bubble));
    trow = trow + offset(1);
    tcol = tcol + offset(2);
     
    % center patch on target location within face and compute intersection
    patch_dims   = RectOfMatrix(patch);
    face_dims    = RectOfMatrix(face);
    t_coords = CenterRectOnPoint(patch_dims, tcol, trow); % left top right bottom, in xy-coordinates!
    t_fit    = ClipRect(face_dims, t_coords);    

    % cut patch to fit into intersection
    start_patch   = (t_fit(1:2) - t_coords(1:2)) + patch_dims(1:2) + 1;
    stop_patch    = (t_fit(3:4) - t_coords(3:4)) + patch_dims(3:4);

    % return values
    bubble_center(bubble,:)  = [trow, tcol];
    bubble_dims(bubble, :)   = [start_patch(2), stop_patch(2), start_patch(1), stop_patch(1)];
    face_coords(bubble, :)   = [t_fit(2) + 1, t_fit(4), t_fit(1) + 1, t_fit(3)];
    end
else
    bubble_center = [];
    bubble_dims   = [];
    face_coords   = [];
end
end
% 
% 
% function [bubble_center, bubble_dims, face_coords] = prepare_alpha(patch, face, n_bubbles);
% % generate random locations for bubbles (center)
% if n_bubbles > 0
% target_center = randperm(numel(face), n_bubbles);
%     
%     for bubble = 1:n_bubbles
%     % location in face where bubble is to be centered
%     [trow, tcol]  = ind2sub(size(face), target_center(bubble));
%     
%     % center patch on target location within face and compute intersection
%     patch_dims   = RectOfMatrix(patch);
%     face_dims    = RectOfMatrix(face);
%     t_coords = CenterRectOnPoint(patch_dims, tcol, trow); % left top right bottom, in xy-coordinates!
%     t_fit    = ClipRect(face_dims, t_coords);    
% 
%     % cut patch to fit into intersection
%     start_patch   = (t_fit(1:2) - t_coords(1:2)) + patch_dims(1:2) + 1;
%     stop_patch    = (t_fit(3:4) - t_coords(3:4)) + patch_dims(3:4);
% 
%     % return values
%     bubble_center(bubble,:)  = [trow, tcol];
%     bubble_dims(bubble, :)   = [start_patch(2), stop_patch(2), start_patch(1), stop_patch(1)];
%     face_coords(bubble, :)   = [t_fit(2) + 1, t_fit(4), t_fit(1) + 1, t_fit(3)];
%     end
% else
%     bubble_center = [];
%     bubble_dims   = [];
%     face_coords   = [];
% end
% end

function [alpha_plane] = add_alphaplane(patch, face, bubble_center, bubble_dims, face_coords);
    alpha_plane = zeros(size(face));
    if ~isempty(bubble_center)
        for bubble = 1:size(bubble_center,1)
          resized_patch = patch(bubble_dims(bubble, 1):bubble_dims(bubble, 2),...
                                bubble_dims(bubble, 3):bubble_dims(bubble, 4));
          tmp = zeros(size(face));
              tmp(face_coords(bubble, 1): face_coords(bubble, 2),...
                  face_coords(bubble, 3):face_coords(bubble, 4)) = resized_patch;    

        alpha_plane = alpha_plane + tmp;
        end
    end
end

function [TastenVector, response_instruction, response_keys] = prepare_responseKeys(response_mapping)
KbName('UnifyKeyNames');
TastenCodes  = KbName({'y','m', 'ESCAPE'}); 
TastenVector = zeros(1,256); TastenVector(TastenCodes) = 1;
if response_mapping == 1
respk_neutral = KbName('y');
respk_emotion = KbName('m');
response_instruction = 'y = Neutral, m = Emotion';
elseif mapping == 2
respk_neutral = KbName('m');
respk_emotion = KbName('y');
response_instruction = 'y = Emotion, m = Neutral';
end
response_keys = [respk_neutral, respk_emotion]; % 1 = neutral, 2 = emotion; for response classification
                                                % response_keys(stim.cond(trial)) == actual_response;
end
