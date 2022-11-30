% dot tracking stimulus
% behavoral experiment with adaptive procedure
% uses 32-bit color images specified by user
% requires psychophysics toolbox routines and statistics toolbox (normrnd function)
% written by Keith Schneider, October 2002, modified from dotrack_behav.m (method of constant stimuli)
% modified 09/18/06 - ported to Psychophysics Toolbox OSX - Mithun Mukherjee
% modified 12/15/14 - changed to an MCS procedure (from a staircase procedure) 


%% This is just a practice program for subject understand the task,total is 10 trials
% 12/29/2014 By modified by Ruyuan Zhang


%%

clear all;close all;clc;
kbindex = -1;KbName('UnifyKeyNames');
% -----------------------
% experimental parameters
% -----------------------

p.total_dots = 8;		% total number of dots

% parameters for MCS procedure
p.attend_dots_cond   = [1,2,3,4,5];
p.trials_per_cond    = [2,2,2,2,2];%trial number corresponding to each level


quitKey              = KbName('Escape');
blueKey              = KbName('LeftArrow');
yellowKey            = KbName('RightArrow');
rkeys                = [blueKey yellowKey]; % keys for "blue face" and "yellow face" responses, respectively (B and Y)


p.attend_dots=[];
for i=1:length(p.attend_dots_cond)
    p.attend_dots=[p.attend_dots repmat(p.attend_dots_cond(i), 1, p.trials_per_cond(i))];
end
p.attend_dots = Shuffle(p.attend_dots);
p.max_trials = length(p.attend_dots);

% % parameters for adaptive procedure
% p.max_trials = 72;		% maximum number of trials -- stop after this many trials even if threshold not found
% p.init_attend_dots = 1;	% initial number of dots to track
% p.min_attend_dots = 1;	% minimum number of dots to track
% p.max_attend_dots = 8;	% maximum number of dots to track (should be <= p.total_dots/2)
% p.correct_inc = 3;		% increase number of dots to track if this many number of correct trials in a row
% p.incorrect_dec = 1;	% decrease number of dots if this many number of incorrect trials in a row
% p.step = 1;				% increment or decrement step
% p.stop_reversals = 8;	% stop after this many reversals (a reversal is a change from increase to decrease or decrease to increase)


p.t_cue = 2;			% time cue is on in seconds (if 0, on all the time)
p.t_move = 4;			% movement time after cue turns off
p.dot_vel = 5;			% velocity of dots in deg/sec

p.move_while_cued = 1;	% should the dots move while the cue is on?

p.circular_field = 1;	% field is circular? (if 0, field is rectangular)

p.p_move_for = 0.4;		% probability that each dot moves straight
p.ang_sd = 0.2;			% sd of angular change profile (the larger this is, the more random the motion appears)

qresponse = 2;			% 1 for clicking mouse on all
						% 2 for highlighting tracked dot (50%) or untracked dot (50%) and key response

qtextoutput = 1;		% output results to text file?


% ------------------------
% set display information
% ------------------------

nScreen = 0;						% which Screen (0 is main Screen)
p.dims.Screen_width = 	40;		    % width of Screen (cm)
p.dims.distance = 		58;			% distance between subject and Screen (cm)

% ------------------------------------
% define Screen and stimuli dimensions
% ------------------------------------

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

p.image_name{1} = 'happy_face.jpg';	% cued
p.image_name{2} = 'sad_face.jpg';	% normal
p.image_name{3} = 'query.jpg';		% query

qthreshbg = 1;	% 1: change colors in image below bgthresh to background field color
				% 2: change colors in image equal to bgthresh to background field color
bgthresh = [0 0 0];	% background threshold or color to make transparent

% -------------------------
% input subject information
% -------------------------

tim=fix(clock);
% p.subject = input('Enter subject''s initials: ', 's');
% if isempty(p.subject)
% 	p.subject = 'test';
% end	

% -------------
% set up trials
% -------------

%p.attend_dots = zeros(p.max_trials,1);
%p.attend_dots(1) = p.init_attend_dots;
p.rtime = zeros(p.max_trials,1);
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
p.correct = zeros(p.max_trials,1);					% 1 if answer is correct, 0 if incorrect	

% -------------
% set up colors
% -------------

col.black = [0 0 0];
col.white = [255 255 255];
col.gray =	[128 128 128];
col.bg =	col.black;	% background color
col.field = col.gray;	% field color
col.done = [0 0 255];

% ---------------
% open the Screen
% ---------------

%%p.res = Screen(nScreen, 'Resolution');
%%p.res = NearestResolution(nScreen, , p.res.height, p.res.hz, 32);
%%w = Screen(nScreen, 'OpenWindow', col.bg, [], p.res);
[w, rect] = Screen('OpenWindow', nScreen, col.bg);
p.res.width = rect(3);
p.res.height = rect(4);

%%p.frame_rate = Screen(w,'FrameRate', []);
frame_duration = Screen('GetFlipInterval', w);
p.frame_rate = 1/frame_duration;

HideCursor;
xc = p.res.width/2;		% horizontal Screen center
yc = p.res.height/2;	% vertical Screen center

% ---------------------------------------------
% set up Screen coordinates in pixel dimensions
% ---------------------------------------------

% number of pixels subtended by one degree of visual angle at fixation
ppd = pi /180 *  p.res.width* p.dims.distance / p.dims.Screen_width;

vel = p.dot_vel * ppd / p.frame_rate;	% velocity in pixels per frame
cvel = ceil(vel);

rdot = round(p.dims.rdot * ppd);												% cue radius (pixels)
min_sep = max(round(p.dims.min_sep * ppd), ceil(2*sqrt(2)*(cvel+1))+2*rdot+3);	% minimum dot separation (pixels)
min_fix = round(p.dims.min_fix * ppd);											% minimum distance allowed from fixation (deg)

if p.circular_field
	max_fix = round(p.dims.max_fix * ppd);	% maximum distaqnce allowed from fixation (pixels)
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

% create image templates
 for i=1:3
%	woff(i) = Screen(w, 'OpenOffScreenWindow', col.field, (rdot+cvel) * [0 0 2 2]);
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


% ----------------
% start experiment
% ----------------

% n_revs = 0;
% step_rising = 0;
% step_falling = 0;
% trials_at_level = 1;
p.final_trial = p.max_trials;

% set up first trial
% draw background field
if p.circular_field
    Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
else
    Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
end
% draw fixation point
Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);


if qresponse == 1
    txt = 'Click the mouse to proceed';
elseif qresponse == 2
    txt = 'Press any key to proceed';
end
Screen('TextSize', w, 15);
wtxt = Screen('TextBounds', w, txt);
txtloc = [xc - wtxt(3)/2, yc + 40];
Screen('DrawText', w, txt, txtloc(1), txtloc(2), col.white);
Screen('Flip', w);
WaitSecs(0.2);

if qresponse == 1
    bdown = 0;
    while(~bdown)
        [x, y, bdown] = GetMouse(w);
    end
    while(bdown)	% wait for mouse button to be released
        [x, y, bdown] = GetMouse(w);
    end
elseif qresponse == 2
    KbWait(kbindex);
end


for t=1:p.max_trials
	
% draw background field
if p.circular_field
 		Screen('FillOval', w, col.field, [xc-max_fix, yc-max_fix, xc+max_fix, yc+max_fix]);
 	else
 		Screen('FillRect', w, col.field, [xc yc xc yc] + [-field_w, -field_h, field_w, field_h]/2);
end
	% draw fixation point
	Screen('FillRect', w, col.white, [xc-2 yc-2 xc+2 yc+2]);
	Screen('FillRect', w, col.black, [xc-1 yc-1 xc+1 yc+1]);
   

% 	if qresponse == 1
% 		txt = 'Click the mouse to proceed';
% 	elseif qresponse == 2
% 		txt = 'Press any key to proceed';
%     end
%     Screen('TextSize', w, 15);
% 	wtxt = Screen('TextBounds', w, txt);
% 	txtloc = [xc - wtxt(3)/2, yc + 40];
%    	Screen('DrawText', w, txt, txtloc(1), txtloc(2), col.white);
     Screen('Flip', w);
 	WaitSecs(1);
% 	
% 	if qresponse == 1
% 		bdown = 0;
% 		while(~bdown)
% 			[x, y, bdown] = GetMouse(w);
% 		end
% 		while(bdown)	% wait for mouse button to be released
% 			[x, y, bdown] = GetMouse(w);
% 		end	
% 	elseif qresponse == 2
% 		KbWait(1);
% 	end


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
	
	rand('state',sum(100*clock));	% initialize random number generator
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
			if ~iflag & i>1
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
	dot_color(1:p.attend_dots(t)) = 2 * ones(1,p.attend_dots(t));
	
	lostframes = 0;
	
	cue_flag = 1;
	
	Screen(w,'WaitBlanking');
	start_t = GetSecs;
	if ~p.move_while_cued
		for i=1:p.total_dots
%		Screen('CopyWindow', woff(dot_color(i)), w, [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1]);		
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
	
%		frames_since_last = Screen(w,'WaitBlanking');
        %Screen('Flip', w); 
%		lostframes = lostframes + frames_since_last;
		
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
			%Screen('CopyWindow', woff(dot_color(i)), w, [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1]);		
            Screen('DrawTexture', w, woff(dot_color(i)), [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1], 0);

        end
        Screen('Flip', w);

		% change direction
		tmp = find(rand(1,p.total_dots) > p.p_move_for);
		if any(tmp)
			%mov_ang(tmp) = mov_ang(tmp) + normrnd(0,p.ang_sd,1,length(tmp));
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
			if sum(dot_state==3) == attend_dots & (xc-x)^2 + (yc-y)^2 < rdot^2	% clicked in center?
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
					if (xpos0(i)-x)^2 + (ypos0(i)-y)^2 <= rdot^2;	% if click located within a dot, select that dot
						
						%Screen('CopyWindow', woff(dot_state(i)), w, [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1]);

                       Screen('DrawTexture', w, woff(dot_state(i)), [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1], 0);

                        dot_state(i) = 4-dot_state(i);

						if dswitch
                            Screen('FillRect', w, field_index, [xc-rdot, yc-rdot, xc+rdot, yc+rdot]);
							Screen('FillRect', w, 255, [xc-2 yc-2 xc+2 yc+2]);
							Screen('FillRect', w, 0, [xc-1 yc-1 xc+1 yc+1]);
							dswitch = 0;
						elseif sum(dot_state==3) == attend_dots
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
		%Screen('CopyWindow', woff(3), w, [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1]);	

        for i=1:p.total_dots
			%Screen('CopyWindow', woff(dot_color(i)), w, [], [xpos0(i), ypos0(i), xpos0(i), ypos0(i)] + (rdot+cvel) * [-1, -1, 1, 1]);		
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
            
            %check if trying to quit out of task immediately
            if (keyCode(quitKey))
                Screen('CloseAll');
                return;
                error('esc:exit','ESC pressed.');
            end
            
			keyresp = keyCode(rkeys);
			if sum(keyresp)  ~= 1	% make sure only one key is pressed
				keyIsDown = 0;
			end
        end
            
		p.answer(t) = find(keyresp);
		p.rtime(t) = secs - rt_start;
		
		p.correct(t) = (p.answer(t) == 1 & p.probe_tracked(t) == 1) | (p.answer(t) == 2 & p.probe_tracked(t) == 0);
		
	end

	p.dot_xpos(t,:) = (xpos0 - xc) / ppd;
	p.dot_ypos(t,:) = (ypos0 - yc) / ppd;
	
	
% 	% adaptive procedure to find level for next trial
% 	if t< p.max_trials
% 		if trials_at_level >= p.correct_inc & prod(p.correct(t-p.correct_inc+1:t)) ...
% 			& p.attend_dots(t) < p.max_attend_dots		% increment level
% 			p.attend_dots(t+1) = p.attend_dots(t) + p.step;
% 			trials_at_level = 1;
% 			step_rising = 1;
% 			if step_falling	% reversal?
% 				step_falling = 0;
% 				n_revs = n_revs + 1;
% 			end
% 		elseif trials_at_level >= p.incorrect_dec & ~sum(p.correct(t-p.incorrect_dec+1:t)) ...
% 			& p.attend_dots(t) > p.min_attend_dots	% decrement level
% 			p.attend_dots(t+1) = p.attend_dots(t) - p.step;
% 			trials_at_level = 1;
% 			step_falling = 1;
% 			if step_rising	% reversal?
% 				step_rising = 0;
% 				n_revs = n_revs + 1;
% 			end
% 		else
% 			p.attend_dots(t+1) = p.attend_dots(t);
% 			trials_at_level = trials_at_level + 1;
% 		end
% 	end
% 
% 	if n_revs == p.stop_reversals
% 		break
% 	end

end

Screen('CloseAll')

% resize arrays to match number of trials completed
p.attend_dots = p.attend_dots(1:t);
p.correct = p.correct(1:t);
p.rtime = p.rtime(1:t);
p.dot_xpos = p.dot_xpos(1:t,:);
p.dot_ypos = p.dot_ypos(1:t,:);
if qresponse == 1
	p.n_attend_clicked = p.n_attend_clicked(1:t);
	p.answer = p.anster(1:t,:);
elseif qresponse == 2
	p.answer = p.answer(1:t);
	p.probed_dot = p.probed_dot(1:t);
	p.probe_tracked = p.probe_tracked(1:t);
end

% sort and score trials by level
p.levels_tested = unique(p.attend_dots);
nlevels = length(p.levels_tested);
p.trials_per_level = zeros(nlevels,1);
p.n_correct_per_level = zeros(nlevels,1);
for i=1:nlevels
	tmp = find(p.attend_dots==p.levels_tested(i));
	p.trials_per_level(i) = length(tmp);
	p.n_correct_per_level(i) = sum(p.correct(tmp));
end
p.mean_p_correct_per_level = p.n_correct_per_level./p.trials_per_level;
p.stderr_p_correct_per_level = 100*sqrt(p.mean_p_correct_per_level.*(1-p.mean_p_correct_per_level)./p.trials_per_level);
p.mean_p_correct_per_level = 100 * p.mean_p_correct_per_level;

% % ----
% % save
% % ----
% 
% file = sprintf('dot-%s%4d%2d%2d%2d%2d%2d', p.subject, tim(1), tim(2), tim(3), tim(4), tim(5), tim(6));
% tmp = find(file == ' ');
% file(tmp) = '0' * ones(size(tmp));
% save(file, 'file', 'p');
% 
% if qtextoutput & qresponse == 2	% output results report to text file
% 	
% 	% compile result codes	
% 	result_letters = ['h','m','c','f']; 	% hit, miss, correct rejection, false alarm
% 	result_nums = zeros(t,1);
% 	tmp = find(p.answer == 1 & p.probe_tracked == 1);
% 	result_nums(tmp) = ones(size(tmp));
% 	tmp = find(p.answer == 2 & p.probe_tracked == 1);
% 	result_nums(tmp) = 2 * ones(size(tmp));
% 	tmp = find(p.answer(1:t) == 2 & p.probe_tracked == 0);
% 	result_nums(tmp) = 3 * ones(size(tmp));
% 	tmp = find(p.answer == 1 & p.probe_tracked == 0);
% 	result_nums(tmp) = 4 * ones(size(tmp));
% 	result_codes = result_letters(result_nums);
% 	
% 	ynletters = ['b','y'];	% letters indicating for expected and actual responses
% 	expected_response = ynletters(2-p.probe_tracked);
% 	actual_response = ynletters(p.answer);
% 	f = fopen([file,'.txt'],'w');
% 	fprintf(f,'TrialNumber CuedDots ExpectedResponse ActualResponse Correct ResultCode\n');
% 	for i=1:t
% 		fprintf(f,'%d %d %s %s %d %s\n',i,p.attend_dots(i),expected_response(i),actual_response(i),p.correct(i),result_codes(i));
% 	end
% 	fclose(f);
% 	
% end

% ------------
% plot results
% ------------

% percent correct per level
figure
errorbar(p.levels_tested, p.mean_p_correct_per_level, p.stderr_p_correct_per_level, p.stderr_p_correct_per_level, 'k-o')
ylabel('Percent correct')	
xlabel('Number of attended dots')
axes('Position',[0.07 0.98 0.15 1],'Visible','off');
%text(0,0,sprintf('%s',file));	% label whole figure

% % levels tested on each trial
% figure
% plot(p.attend_dots,'b-x')
% xlabel('Trial number')
% ylabel('Number of dots to track')
% axes('Position',[0.07 0.98 0.15 1],'Visible','off');
% text(0,0,sprintf('%s',file));	% label whole figure