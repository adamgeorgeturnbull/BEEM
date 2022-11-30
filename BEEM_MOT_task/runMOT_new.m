function runMOT_new(subID,session_number,PoF,num_dots_track)
% Inputs:
% subID - subject's ID. Must be a number.
%
% session_number - session number for the given subject.
%
% PoF - an indicator variable for whether the current run is a practice
% run (PoF = 0) or the actual experiment run (PoF = 1). No data will
% be saved from a practice run.
%
% num_dots_track - number of faces that subject will track during each
% trial. Must be an integer smaller than the total number of faces which
% move during each trial.

%% Preliminary Steps

%--------------------------------------------------------------------------
% cd to correct Desktop folder
comp_type = computer;
if strcmp(comp_type(1:3),'MAC') == 1
    desktop_path = fullfile('~','Desktop');
elseif strcmp(comp_type(1:5),'PCWIN') == 1
    desktop_path = winqueryreg('HKEY_CURRENT_USER', 'Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders', 'Desktop'); 
end

cd(fullfile(desktop_path,'MOT_task_materials'))

%--------------------------------------------------------------------------
% Check if this is a training session or an actual session. No data is
% collected from practice sessions, while data from actual sessions will be
% saved.
if PoF == 0
    test_yn = 1;
elseif PoF == 1
    test_yn = 0;
end

% Create folder for subject.
% results_folder = fullfile(pwd,'results');
results_folder = fullfile(pwd,'tDCS_study_MOT_results');

% Get time of day when task begins.
full_clock = clock;
hours = full_clock(4);
minutes = full_clock(5);
seconds = round(full_clock(6));

if hours < 10
    hr_string = ['0' num2str(hours)];
else
    hr_string = num2str(hours);
end

if minutes < 10
    min_string = ['0' num2str(minutes)];
else
    min_string = num2str(minutes);
end

if seconds < 10
    sec_string = ['0' num2str(seconds)];
else
    sec_string = num2str(seconds);
end

task_start_time = [hr_string '-' min_string '-' sec_string];

% Create dataset name that includes subject ID, session number, date, and
% time of day when task began.

dataset_name = ['ID' num2str(subID) '_Session' num2str(session_number) '_' date '_' task_start_time];

if test_yn == 0
    sub_folder = fullfile(results_folder,['ID' num2str(subID)]);
    if ~exist(sub_folder,'dir')
        mkdir(sub_folder)
    end
%     sub_session_folder = fullfile(sub_folder,['Session' num2str(session_number) '_' date '_' task_start_time]);
%     if ~exist(sub_session_folder,'dir')
%         mkdir(sub_session_folder)
%     else
%         error('A folder for subject %i session %i already exists. Delete that folder and re-run runMOT_new if you want to overwrite that data.',subID,session_number)
%     end

    dataset_path = fullfile(sub_folder,dataset_name);
end


%--------------------------------------------------------------------------
% Setup colors
col.black = [0 0 0];
col.white = [255 255 255];
col.gray =	[128 128 128];
col.bg =	[255 255 255];	% background color
col.field = [153 205 255];	% field color
col.done = [0 0 255];

%--------------------------------------------------------------------------
% Setup screen and keyboard

% Select screen.
nScreen = max(Screen('Screens'));	

% Screen setup
clear screen
Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference', 'VisualDebugLevel', 0);
Screen('Preference', 'SkipSyncTests', 1);
[w, rect] = Screen('OpenWindow', nScreen, col.bg, [], [], 2);
p.res.width = rect(3);
p.res.height = rect(4);

frame_duration = Screen('GetFlipInterval', w);
p.frame_rate = 1/frame_duration; % Get monitor's frame rate.

HideCursor;
xc = p.res.width/2;		% horizontal Screen center
yc = p.res.height/2;	% vertical Screen center

% Keyboard Setup
clear PsychHID; % Force new enumeration of devices.
clear KbCheck; % Clear persistent cache of keyboard devices. 
KbName('UnifyKeyNames');    % keynames for all types of keyboard
kbindex = -3;

%--------------------------------------------------------------------------
% Setup experiment parameters
p.total_dots = 8;		% total number of dots

% parameters for MCS procedure
p.attend_dots_cond   = [num_dots_track,num_dots_track,num_dots_track,num_dots_track];
p.trials_per_cond    = [50,200,200,200,200];%trial number corresponding to each level


quitKey              = KbName('Escape');
blueKey              = KbName('8');
yellowKey            = KbName('9');
rkeys                = [blueKey yellowKey]; % keys for "blue face" and "yellow face" responses, respectively (B and Y)

p.attend_dots=[];
for i=1:length(p.attend_dots_cond)
    p.attend_dots=[p.attend_dots repmat(p.attend_dots_cond(i), 1, p.trials_per_cond(i))];
end
p.attend_dots = Shuffle(p.attend_dots);%repmat(num_dots_track,[1 length(p.attend_dots_cond)]);
p.max_trials = length(p.attend_dots);

% Timing and movement parameters
p.t_cue = 2;			% time cue is on in seconds (if 0, on all the time)
p.t_move = 4;			% movement time after cue turns off
p.dot_vel = 5;			% velocity of dots in deg/sec

p.move_while_cued = 1;	% should the dots move while the cue is on?

p.circular_field = 1;	% field is circular? (if 0, field is rectangular)

p.p_move_for = 0.4;		% probability that each dot moves straight
p.ang_sd = 0.2;			% sd of angular change profile (the larger this is, the more random the motion appears)

qresponse = 2;			% 1 for clicking mouse on all, 2 for highlighting tracked dot (50%) or untracked dot (50%) and key response

%--------------------------------------------------------------------------
% define Screen size
p.dims.Screen_width = 	40;		    % width of Screen (cm)
p.dims.distance = 		58;			% distance between subject and Screen (cm)

% Provide stimulus dimensions (in degrees)
p.dims.rdot = 			0.4;		% dot radius (deg)
p.dims.min_sep =		1.5;		% minimum distance allowed between dots (deg)
p.dims.min_fix =		3;			% minimum distance allowed from fixation (deg)
if p.circular_field
    p.dims.max_fix =	10;			% maximum distance allowed from fixation (deg)
else
    p.dims.field_h =	20;			% height of field (deg)
    p.dims.field_w =	20;			% width of field (deg)
end
p.dims.min_edge =		0;			% minimum distance from edge (deg)

%--------------------------------------------------------------------------
% Convert screen coordinates to pixels
% number of pixels subtended by one degree of visual angle at fixation
ppd = pi /180 *  p.res.width* p.dims.distance / p.dims.Screen_width;

vel = p.dot_vel * ppd / p.frame_rate;	% velocity in pixels per frame
cvel = ceil(vel);

rdot = round(p.dims.rdot * ppd);												% cue radius (pixels)
min_sep = max(round(p.dims.min_sep * ppd), ceil(2*sqrt(2)*(cvel+1))+2*rdot+3);	% minimum dot separation (pixels)
min_fix = round(p.dims.min_fix * ppd);											% minimum distance allowed from fixation (deg)

if p.circular_field
    max_fix = round(p.dims.max_fix * ppd);	% maximum distance allowed from fixation (pixels)
    min_edge = max(round(p.dims.min_edge * ppd), ceil(2*sqrt(2)*(cvel+1))+rdot+4);		% minimum distance from edge (pixels)
else
    min_edge = max(round(p.dims.min_edge * ppd), rdot + cvel + 2);		% minimum distance from edge (pixels)
    field_h = round(p.dims.field_h * ppd);	% field height (pixels)
    field_w = round(p.dims.field_w * ppd);	% field width (pixels)
end


if p.circular_field
    x0 = min_edge;
    x1 = p.res.width-min_edge;
    y0 = min_edge;
    y1 = p.res.height-min_edge;
else
    x1 = xc+field_w/2 - min_edge;
    x0 = xc-field_w/2 + min_edge;
    y1 = yc+field_h/2 - min_edge;
    y0 = yc-field_h/2 + min_edge;
end
dstep = pi/50;	% step per frame by which to change direction angle towards confined quadrant

bounds = [x0 x1 y0 y1];
boundn = ones(1,p.total_dots);

%--------------------------------------------------------------------------
% Load images and create image templates.
p.image_name{1} = 'happy_face.png';	% cued
p.image_name{2} = 'sad_face.png';	% normal
p.image_name{3} = 'query.png';		% query
image_correct = imread('Correct.png');
image_incorrect = imread('Incorrect.png');

qthreshbg = 1;	% 1: change colors in image below bgthresh to background field color
% 2: change colors in image equal to bgthresh to background field color
bgthresh = [0 0 0];	% background threshold or color to make transparent

for i=1:numel(p.image_name)
    tmp = imread(p.image_name{i});
    if qthreshbg	% make parts of image transparent?
        if qthreshbg == 1
            tothresh = find(tmp(:,:,1) <= bgthresh(1) & tmp(:,:,2) <= bgthresh(2) & tmp(:,:,3) <= bgthresh(3));
        else
            tothresh = find(tmp(:,:,1) == bgthresh(1) & tmp(:,:,2) == bgthresh(2) & tmp(:,:,3) == bgthresh(3));
        end
        for j=1:3
            tmp1 = tmp(:,:,j);
            tmp1(tothresh) = col.field(j)*ones(size(tothresh));
            tmp(:,:,j) = tmp1;
        end
    end
    %	Screen(woff(i),'PutImage', tmp, cvel + 2*rdot * [0 0 1 1]);
    woff(i) = Screen('MakeTexture', w, tmp);
end

%--------------------------------------------------------------------------
% Setup behavioral datasets
p.rtime = zeros(p.max_trials,1);
p.trial_start = zeros(p.max_trials,1);
if qresponse == 1
    p.answer = zeros(p.max_trials, p.total_dots);
    p.n_attend_clicked = zeros(p.max_trials,1);
else
    p.answer = zeros(p.max_trials, 1); % 1 for yes 2 for no
    p.probe_tracked = floor(2*rand(p.max_trials,1));
    p.probed_dot = zeros(p.max_trials,1);				% the number of the probed dot
end
p.dot_xpos = zeros(p.max_trials, p.total_dots);
p.dot_ypos = zeros(p.max_trials, p.total_dots);
p.correct = NaN*ones(p.max_trials,1);					% 1 if answer is correct, 0 if incorrect

%% Begin Trials
%--------------------------------------------------------------------------
% Setup first screen
p.final_trial = p.max_trials;

% draw background field
if p.circular_field
    Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
else
    Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
end

% draw fixation point
Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);

% First screen text
if qresponse == 1
    txt = 'Click the mouse to proceed';
elseif qresponse == 2
    txt = 'Please wait...';
end

% Display screen with text and fixation point over background field.
Screen('TextSize', w, 40);
DrawFormattedText(w,txt,'center',0.40*p.res.height,col.white)
Screen('Flip', w);
WaitSecs(0.2);

% Wait for administrator to click mouse to begin trials.
if qresponse == 1
    bdown = 0;
    while(~bdown)
        [x, y, bdown] = GetMouse(w);
    end
    while(bdown)	% wait for mouse button to be released
        [x, y, bdown] = GetMouse(w);
    end
elseif qresponse == 2
    % Wait for the interviewer pressing the mouse
    while 1
        [clicks,x,y,whichButton] = GetClicks;
        if clicks>0
            break
        end
    end
end

%--------------------------------------------------------------------------
% Run trials
trials_total = 0;
tic; % record the start time of task

for t=1:p.max_trials
    trial_start = toc; % record the start time for each trial
    trials_total = trials_total+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Draw background field and fixation point, hold for 1 second.
    % draw background field
    if p.circular_field
        Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
    else
        Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
    end
    
    % draw fixation point
    Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
    Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);
    Screen('Flip', w);
    WaitSecs(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % draw background field
    if p.circular_field
        Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
    else
        Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
    end
    % draw fixation point
    Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
    Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);
    
    % --------------------------------------
    % chose initial positions and velocities
    % --------------------------------------
    
    xpos0 = zeros(1,p.total_dots);
    ypos0 = xpos0;
    % choose initial positions but make sure they are separated at least by min_sep pixels
    
    rng(sum(100*clock));	% initialize random number generator
    for i = 1:p.total_dots
        iflag = 1;
        while(iflag)
            iflag = 0;
            if p.circular_field
                xpos0(i) = rand(1) * 2 * (max_fix-min_edge) + min_edge + xc - max_fix;
                ypos0(i) = rand(1) * 2 * (max_fix-min_edge) + min_edge + yc - max_fix;
            else
                xpos0(i) = rand(1) * (field_w-2*min_edge) + min_edge + xc-field_w/2;
                ypos0(i) = rand(1) * (field_h-2*min_edge) + min_edge + yc-field_h/2;
            end
            r2 = (xpos0(i)-xc)^2 + (ypos0(i)-yc)^2;
            if r2 < min_fix^2
                iflag = 1;
            elseif p.circular_field
                if r2 > (max_fix-min_edge)^2
                    iflag = 1;
                end
            end
            if ~iflag && i>1
                for j = 1:i-1
                    if (xpos0(i)-xpos0(j))^2 + (ypos0(i)-ypos0(j))^2 < min_sep^2
                        iflag = 1;
                        break
                    end
                end
            end
        end
    end
    mov_ang = rand(1,p.total_dots) * 2 * pi;
    
    dot_color = ones(1,p.total_dots);
    dot_color(1:p.attend_dots(t)) = 2 * ones(1,p.attend_dots(t));%p.attend_dots(t)
    
    lostframes = 0;
    
    cue_flag = 1;
    
    Screen('WaitBlanking',w);
    start_t = GetSecs;
    if ~p.move_while_cued
        for i=1:p.total_dots
            Screen('DrawTexture', w, woff(dot_color(i)), [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1], 0);
        end
        
        WaitSecs(p.t_cue);
        cue_flag = 0;
        dot_color = ones(1,p.total_dots);
        Screen('Flip', w);
    end
    
    
    while(GetSecs - start_t < p.t_move + p.t_cue)
        if cue_flag
            if GetSecs - start_t >= p.t_cue
                cue_flag = 0;
                dot_color = ones(1,p.total_dots);
            end
        end
        
        % draw background field
        if p.circular_field
            Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
        else
            Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
        end
        % draw fixation point
        Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
        Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);
        
        for i=1:p.total_dots
            Screen('DrawTexture', w, woff(dot_color(i)), [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1], 0);
        end
        Screen('Flip', w);
        
        % change direction
        tmp = find(rand(1,p.total_dots) > p.p_move_for);
        if any(tmp)
            r = randn(1, length(tmp)) .* p.ang_sd;
            mov_ang(tmp) = mov_ang(tmp) + r;
        end
        
        xpos = xpos0 + cos(mov_ang) * vel;	% predicted position change
        ypos = ypos0 - sin(mov_ang) * vel;
        
        % reflect off horizontal sides
        tmp = find(xpos < bounds(boundn, 1)' | xpos > bounds(boundn, 2)');
        if any(tmp)
            mov_ang(tmp) = pi - mov_ang(tmp);
        end
        %reflect off vertical sides
        tmp = find(ypos < bounds(boundn, 3)' | ypos > bounds(boundn, 4)');
        if any(tmp)
            mov_ang(tmp) = 2*pi - mov_ang(tmp);
        end
        % reflect off center or outer circle
        r2 = (xpos-xc).^2 + (ypos-yc).^2;
        if p.circular_field
            tmp = find(r2 < min_fix^2 | r2 > (max_fix-min_edge)^2);
        else
            tmp = find(r2 < min_fix^2);
        end
        if any(tmp)
            mov_ang(tmp) = 2 * atan2(-(ypos0(tmp)-yc),xpos0(tmp)-xc) - mov_ang(tmp) - pi;
        end
        % reflect dots off each other
        for i=1:p.total_dots-1
            for j = i+1:p.total_dots
                if (xpos(i)-xpos(j))^2 + (ypos(i)-ypos(j))^2 < min_sep^2
                    mov_ang([i j]) = mov_ang([j i]);
                end
            end
        end
        xpos0 = xpos0 + cos(mov_ang) * vel;	% actual position change
        ypos0 = ypos0 - sin(mov_ang) * vel;
    end
    
    xpos0 = xpos0 - cos(mov_ang) * vel;	% change position back at end of trial
    ypos0 = ypos0 + sin(mov_ang) * vel;
    
    if qresponse == 1	% click mouse on all of tracked dots
        
        setmouse(xc,yc,w);
        ShowCursor(0)
        
        done = 0;
        dot_state = ones(1,p.total_dots);
        
        % click on the cued dots; double-click on center when done
        
        bdown = 0;
        dswitch = 0;
        
        rt_start = GetSecs;
        while(~done)
            while(~bdown)
                [x, y, bdown] = GetMouse(w);
            end
            if sum(dot_state==3) == p.attend_dots && (xc-x)^2 + (yc-y)^2 < rdot^2	% clicked in center?
                done = 1;	% done if click in center
            else
                % draw background field
                if p.circular_field
                    Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
                else
                    Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
                end
                % draw fixation point
                Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
                Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);
                for i=1:p.total_dots
                    if (xpos0(i)-x)^2 + (ypos0(i)-y)^2 <= rdot^2	% if click located within a dot, select that dot
                        
                        %Screen('CopyWindow', woff(dot_state(i)), w, [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1]);
                        
                        Screen('DrawTexture', w, woff(dot_state(i)), [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1], 0);
                        
                        dot_state(i) = 4-dot_state(i);
                        
                        if dswitch
                            Screen('FillRect', w, field_index, [xc-rdot, yc-rdot, xc+rdot, yc+rdot]);
                            Screen('FillRect', w, 255, [xc-2 yc-2 xc+2 yc+2]);
                            Screen('FillRect', w, 0, [xc-1 yc-1 xc+1 yc+1]);
                            dswitch = 0;
                        elseif sum(dot_state==3) == p.attend_dots
                            Screen('FillOval', w,  col.done, [xc-rdot, yc-rdot, xc+rdot, yc+rdot]);
                            Screen('DrawText', w, 'Done?',xc-17,yc+4,col.white);
                            dswitch = 1;
                        end
                        Screen('Flip', w);
                        break
                    end
                end
            end
            while(bdown)	% wait for mouse button to be released
                [x, y, bdown] = GetMouse(w);
            end
        end
        p.rtime(t) = GetSecs - rt_start;
        HideCursor
        
        p.answer(t,:) = dot_state == 3;
        p.n_attend_clicked(t) = sum(p.answer(t,1:p.attend_dots(t)));
        p.correct(t) = p.n_attend_clicked(t) == p.attend_dots(t);
        
    elseif qresponse == 2
        
        if p.probe_tracked(t)
            p.probed_dot(t) = floor(p.attend_dots(t)*rand)+1;
        else
            p.probed_dot(t) = floor((p.total_dots-p.attend_dots(t))*rand)+p.attend_dots(t)+1;
        end
        
        % draw background field
        if p.circular_field
            Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
        else
            Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
        end
        % draw fixation point
        Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
        Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);
        
        for i=1:p.total_dots
            Screen('DrawTexture', w, woff(dot_color(i)), [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1], 0);   
        end
        
        i=p.probed_dot(t);
        Screen('DrawTexture', w, woff(3), [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1], 0);
        Screen('Flip', w);
        rt_start = GetSecs;
        
        % get response
        keyIsDown = 0;
        while(~keyIsDown)
            [keyIsDown,secs,keyCode] = KbCheck(kbindex);
            
            %check if administrator wants to quit out of task.
            if (keyCode(quitKey))                
                break;
            elseif keyIsDown && keyCode(quitKey) == 0
                keyresp = keyCode(rkeys);
                if sum(keyresp)  ~= 1	% make sure only one key is pressed
                    keyIsDown = 0;
                else
                    p.answer(t) = find(keyresp);
                    p.rtime(t) = secs - rt_start;
                    p.trial_start(t) = trial_start;
                    p.correct(t) = (p.answer(t) == 1 & p.probe_tracked(t) == 1) | (p.answer(t) == 2 & p.probe_tracked(t) == 0);
                    
                    p.dot_xpos(t,:) = (xpos0 - xc) / ppd;
                    p.dot_ypos(t,:) = (ypos0 - yc) / ppd;

                    % show feedback
                    display_correct = Screen('MakeTexture', w, image_correct);
                    display_incorrect = Screen('MakeTexture', w, image_incorrect);
                    feedback_img_size = 0.05*min(p.res.width,p.res.height);

                    if p.correct(t)
                        Screen('DrawTexture', w, display_correct, [], [ (xc-feedback_img_size) (yc-feedback_img_size) (xc+feedback_img_size) (yc+feedback_img_size)]);
                        feedback = Screen('Flip', w);
                    else
                        Screen('DrawTexture', w, display_incorrect, [], [ (xc-feedback_img_size) (yc-feedback_img_size) (xc+feedback_img_size) (yc+feedback_img_size)]);
                        feedback = Screen('Flip', w);
                    end
                    WaitSecs(0.5);
                end
            end
        end
    end
        
    % if administrator presses 'Escape' key, end task and save data.
    if (keyCode(quitKey))
        p.trial_start(t) = trial_start;
        break;
    end
end
%% Save behavior data
if trials_total > 1
    trials_total = trials_total - 1;
end

rt_data = p.rtime(1:trials_total);
trial_start_times = p.trial_start(1:trials_total);
response = p.answer(1:trials_total);
probe_tracked_yn = p.probe_tracked(1:trials_total);
probed_dot_label = p.probed_dot(1:trials_total);
correct_yn = p.correct(1:trials_total);
trials = (1:trials_total)';
subID_var = subID*ones(trials_total,1);
num_faces_track = num_dots_track*ones(trials_total,1);

combined_data = horzcat(subID_var,trials,num_faces_track,probed_dot_label,probe_tracked_yn,response,correct_yn,rt_data,trial_start_times);

headers = {'Subject ID','Trial Number','Number of Targets','Probed Object Label','Probed Track (y/n)','Response','Correct (y/n)','Reaction Time','Trial Start Time'};
MOT_behav_data = [headers; num2cell(combined_data)];

if test_yn == 0
%     save(fullfile(sub_session_folder,'MOT_behav_data'),'MOT_behav_data')
    save(dataset_path,'MOT_behav_data')
end

sca;
return;
end