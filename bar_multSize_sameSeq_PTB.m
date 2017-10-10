% bar_multSize_sameSeq_PTB.m
%
% couldn't take the MGL bullshit, doing this in PTB
%
% TCS 6/10/2017
%
% subj: subj name (for ET, only first 3 chars used)
% run:  run number
% seq:  leave out by default; otherwise a n_sweep entry matrix with first
% col for bar width (deg) and second col for bar direction (1/2 = down/up,
%     3/4 = left/right) - matches wayne/masih's [IN FUTURE: 0's in either
%     entry will correspond to blanks]
% ALL COORDS IN CARTESIAN COORD FRAME!!!! + is UP on screen
%
% TODO: get correct dot density, average speed
%
% 1.2 s TR: 304 TRs
% 0.8 s TR: 456 TRs
%
% TCS 7/18/2017 - added eyetracking, added support for behavioral testing



function bar_multSize_sameSeq_PTB(subj,run,seq)

p.subj = subj;
p.run = run;

p.filename = sprintf('./data/%s_r%02.f_RF_bar_multSize_%s.mat',p.subj,p.run,datestr(now,30));
if ~exist('./data','dir')
    mkdir('./data');
end

% seed random number generator (for trial response order, etc)
ctime = cputime*1000;
rng(ctime);

% if demo, do fast stimuli and fast wait at beginning
p.scanner = 1;

p.do_et = 1;

if p.do_et == 1
    p.eyedatafile = sprintf('%s_RF%02.f',p.subj(1:min(length(p.subj),3)),p.run);
end


if nargin < 3
    % made these bigger for scanner - if we're doign 12 steps, to fully
    % sample space at small bar, need ~2.5 deg width; I think mackey et
    % al did sub-sampling of space on smallest bar, though
%  p.seq=[2     2
%         1     3
%         3     2
%         2     4
%         1     2
%         2     3        
%         1     1
%         3     3        
%         3     1
%         1     4                
%         3     4
%         2     1] .* [2.5 1];

% ABOVE REQUIRES 2016b or newer!!!!
    p.seq = [2     2
        1     3
        3     2
        2     4
        1     2
        2     3        
        1     1
        3     3        
        3     1
        1     4                
        3     4
        2     1];
    p.seq(:,1) = p.seq(:,1)*2.5;

else
    p.seq = seq;
end


% ----- bar stimulus info -----

if p.run>1 && exist(sprintf('./data/%s_r%02.f_lastCoh.mat',p.subj,p.run-1),'file')
    load(sprintf('./data/%s_r%02.f_lastCoh.mat',p.subj,p.run-1),'init_coh');
    p.coh_init = init_coh;
    clear init_coh;
else
    p.coh_init    = 0.32; % TODO: load from file....
end

p.n_steps     = 12; % to match previous
p.step_dur    = 2.4;% 2.5; % sec
p.n_segments  = 3; % for now
p.dot_speed   = 1.6; % deg/s
p.dot_life    = 3*2;   % frames
p.coh_sample  = 0.5; % coherence of sample stimulus (outer bars)
p.coh_step    = 0.075; % how much to step up/down in staircase
p.dot_color   = [1 1 1]*255;
p.segment_gap = 0.25; % deg, gap between segments - determines SIZE of rect, will be centered on a known point, so will be gap/2 at top/bottom/left/right too
p.bar_extent  = 1;   % % of square the bar subtends/100 (NOTE: gap will be drawn at top/bottom or left/right of entire bar, too)
p.dot_density = 10; % per deg^2 (seems high...but from biorxiv... - 124 dots per 4 deg^2
p.dot_size_deg= 0.075; % deg

% assign segments (1..n_segments, left->right, top->bottom)
if p.n_segments == 3 % note: only option right now....could get fancier? distractors?
    p.target_segments = [1 3]; % seg 1 is top/left
    p.sample_segments = [2];
end

% ----- display info (unrelated to bar) -------
p.bg_color  = [0.5 0.5 0.5]*255;
p.fb_colors = [0 0.65 0; 0.8 0 0]*255; % correct; incorrect
p.fix_color = [0.75 0.75 0.75] * 255;
p.fix_size  =  0.25; % radius, dva
p.fix_width =  2; % line_width, pixels...
p.fix_aperture_size = 0.75; % radius, dva

% potentially a square aperture? or aperture that looks like bore insert?

if p.scanner == 0
    p.resolution = [1280 1024]; % desired resolution, compare to the 'actual' resolution below
    p.refresh_rate = 60;
    p.screen_height = 30; % cm, in the experiment room
    p.viewing_distance = 56; % cm, in the experiment room (inside lab)
else
    p.resolution = [1280 1024]; % scanner: FOR EYETRACKING!!!
    p.refresh_rate = 120;
    p.screen_height = 36; % cm prior to 9/14/2017 - was 33 cm, adjusted, now 36 cm tall
    p.viewing_distance = 63 + 9.5; % cm
end

p.screen_width = p.screen_height * p.resolution(1)/p.resolution(2); %
p.screen_height_deg = 2*atan2d(p.screen_height/2,p.viewing_distance);
p.screen_width_deg  = 2*atan2d(p.screen_width/2, p.viewing_distance);
p.ppd = p.resolution(2)/p.screen_height_deg;  % used to convert rects, positions later on

p.scr_center = p.resolution/2;  % could do offset centers, etc?

fix_aperture_rect = CenterRectOnPoint(p.ppd * [0 0 2 2] * p.fix_aperture_size,p.scr_center(1),p.scr_center(2));

% ----- compute bar properties along its axis -------
% aperture along long axis of bar for each segment - width determined by
% seq
p.bar_segment_size_deg = (p.bar_extent * p.screen_height_deg/p.n_segments) - p.segment_gap;
p.bar_segment_centers_deg = (((1:p.n_segments)-0.5)-p.n_segments/2) * (p.bar_extent * p.screen_height_deg / p.n_segments); % where the segments are centered




% ----- timing info (unrelated to bar) ------
if p.scanner == 1
    p.start_wait = 9.6; % s - time after trigger
    p.end_wait   = 9.6; % s - time after last bar sweep (2.4 s x 4)
else
    p.start_wait = 1;
    p.end_wait = 1;
end


% ----- things to keep track of --------------

p.bar_pos   = nan(p.n_steps * size(p.seq,1),2); % x, y coord of middle of bar (dva)
p.resp      = nan(p.n_steps * size(p.seq,1),1); % response key (left/right or up/down)
p.correct   = nan(p.n_steps * size(p.seq,1),1); % was response correct?
p.rt        = nan(p.n_steps * size(p.seq,1),1); % when was response relative to beginning of stim?
p.sample_coherence = p.coh_sample*ones(p.n_steps * size(p.seq,1),1); % what was coherence on that stim?
p.target_coherence = nan(p.n_steps * size(p.seq,1),1); % what was coherence on that stim?
p.target_coherence(1) = p.coh_init; % middle bar
p.segment_dirs = nan(p.n_steps * size(p.seq,1),p.n_segments); % direction of each dot patch?

p.segment_dirs(:,[1 2]) = round(rand(size(p.segment_dirs,1),2))*2-1;
p.segment_dirs(:,3) = -1 * p.segment_dirs(:,1);

% staircase
p.up_thresh = 3;   % 3 correct,   make it harder
p.down_thresh = 1; % 1 incorrect, make it easier


% what the correct response is:
p.corr_resp = (p.segment_dirs(:,2) == p.segment_dirs(:,3)) + 1;
% TCS quickf ix
%p.corr_resp(p.seq(:,2)==1|p.seq(:,2)==2) = mod(p.corr_resp(p.seq(:,2)==1|p.seq(:,2)==2),2)+1;


% ----- timing to keep track of --------------
p.step_start  = nan(p.n_steps * size(p.seq,1),1);
p.sweep_start = nan(size(p.seq,1),1); % redundant, but maybe convenient?
p.sweep_end   = nan(size(p.seq,1),1);


% ----- button info ------
if ismac == 1
    p.esc_key = KbName('escape'); % press this key to abort
else
    p.esc_key = KbName('esc');
end
p.start_key = [KbName('5%') KbName('5')];  % should be lower-case t at prisma? (or %5, which is top-row 5, or 5, which is numpad 5)
p.resp_keys = [KbName('1!'),KbName('2@')]; % 1 = left/up, 2 = right/down






% ----- init screen (geometry; optics) ------
% sync
Screen('Preference', 'SkipSyncTests', 0);
if p.scanner == 1
    [w,rect] = Screen('OpenWindow', 2,p.bg_color);
else    
    [w,rect] = Screen('OpenWindow', max(Screen('Screens')),p.bg_color);
end

% refresh rate
ifi = Screen('GetFlipInterval', w);
p.actual_refresh  = 1 / ifi;



Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


p.resolution_actual = [rect(3) rect(4)]; % pixels

HideCursor;

% --------- eyetracking ----------- %
if p.do_et == 1
    
    if p.scanner == 1
        Eyelink('SetAddress','192.168.1.5')
    end
    
    el=EyelinkInitDefaults(w);
    
    el.backgroundcolour=p.bg_color(1);  % TODO: fix this?
    el.calibrationtargetcolour=p.fix_color(1);    
    %el.calibrationtargetwidth=1;
    el.msgfontcolour=p.fix_color(1);
    p.foregroundcolour=p.fix_color(1);
    
    EyelinkUpdateDefaults(el);

    
    Eyelink('Initialize','PsychEyelinkDispatchCallback') % initialises the eyetracker
    
    % SCANNER: right eye!!!!!! TODO: script this, record output...
    Eyelink('command','calibration_type=HV5'); % updating number of callibration dots
    s=Eyelink('command','link_sample_data=LEFT,RIGHT,GAZE,AREA');% (,GAZERES,HREF,PUPIL,STATUS,INPUT');
    s=Eyelink('command', 'sample_rate=1000');
    s=Eyelink('command','screen_pixel_coords=%ld %ld %ld %ld', 0, 0, rect(3)-1,rect(4)-1);
    s=Eyelink('message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, rect(3)-1,rect(4)-1);
    s=Eyelink('command','file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON');
    s=Eyelink('command','file_sample_data = LEFT,RIGHT,GAZE,AREA,GAZERES,STATUS');
    
    
    % make sure that we get gaze data from the Eyelink
    
    
    
    %------ calibrate the eye tracker --------
    EyelinkDoTrackerSetup(el);
    if s~=0
        error('link_sample_data error, status: ',s)
    end
    Eyelink('openfile',p.eyedatafile);

end



% ------ beginning of run: fixation ------


Screen('DrawLines',w, [-1 1 0 0; 0 0 -1 1]*p.ppd*p.fix_size,p.fix_width,0.8*p.fix_color,p.scr_center);
Screen('Flip',w);
% 
% Screen('Preference','TextRenderer',1);
% Screen('DrawLines',w, [-1 1 0 0; 0 0 -1 1]*p.ppd*p.fix_size,p.fix_width,p.fix_color,p.scr_center);
% txt = 'test';%sprintf('Acc: %.02f%%',100*nansum(p.correct)/p.num_trials);
% Screen('TextSize', w, 30);
% DrawFormattedText(w,txt,'center',p.scr_center(2)-3*p.ppd,p.fix_color);
% %normBoundsRect = Screen('TextBounds', w, txt);
% %txtloc = [p.xc - normBoundsRect(3)/2, 2 * p.ppd + p.yc - normBoundsRect(4)/2];
% %Screen('DrawText', w, txt,txtloc(1),txtloc(2), 255 );
% Screen('Flip',w);


% ------ wait for trigger, then brighten fix -------

resp = 0;
fprintf('waiting for start\n');
while resp == 0
    [resp, ~] = checkForResp(p.start_key, p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during pre-trigger wait time\n');
        ShowCursor;
        return;
    end
end
clear resp;

Screen('DrawLines',w, [-1 1 0 0; 0 0 -1 1]*p.ppd*p.fix_size,p.fix_width,p.fix_color,p.scr_center);
Screen('Flip',w);

p.expt_start = GetSecs;

if p.do_et == 1
    Eyelink('Message','xDAT %i', 101);
    Eyelink('StartRecording'); % make 1 big edf file (save time)
end


% ------ initial wait time (check for esc) ---------

resp = 0;
while (GetSecs-p.expt_start) < p.start_wait
    [resp, ~] = checkForResp(p.start_key, p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during post-trigger wait time\n');
        ShowCursor;
        return;
    end
end
clear resp;


% ------ loop over sweeps ---------------------------

trial_counter = 1;

up_cnt = 0; down_cnt = 0;

for ss = 1:size(p.seq,1)
    
    if p.do_et == 1
        
        Eyelink('Message','xDAT %i',ss);
        
        Eyelink('command', 'record_status_message "SWEEP %d of %d"', ss, size(p.seq,1));
        
    end

    
    
    % --------- init sweep ----------
    % ensure bar reaches out to, and only out to, full extent of screen used
    this_bar_pos = linspace(-0.5*p.bar_extent*p.screen_height_deg + p.seq(ss,1)/2, 0.5*p.bar_extent*p.screen_height_deg - p.seq(ss,1)/2, p.n_steps);
    
    
    % if p.seq(ss,2)==1 (down) or 3 (left), flip
    if p.seq(ss,2)==1 || p.seq(ss,2)==3
        this_bar_pos = this_bar_pos(end:-1:1);
    end
    
    % uniform dot grid from which we sample our init dots
    
    % first coord: always along long axis of bar, second coord, parallel to
    % direction of motion (like horiz bar)
    
    segment_area = p.bar_segment_size_deg * p.seq(ss,1);
    ndots = round(p.dot_density * segment_area);
    
    % holds coordinates of dots which are updated on each frame
    bar_dot_coords = cell(p.n_segments,1);
    bar_dot_age = cell(p.n_segments,1);   % current age
    bar_dot_dir = cell(p.n_segments,1);   % direction of each dot (only updated at/after death)
    
    % generate too many dots, to make sure we have enough; if we sample
    % ndots from the correct region, we'll achieve average density over the
    % segment as targeted (p.dot_density, in dots/deg^2)
    tmpndots = round(2*p.dot_density * max(p.bar_segment_size_deg,p.seq(ss,1)).^2);
    for bb = 1:p.n_segments
        coords_tmp = (rand(tmpndots,2)-0.5)*max(p.bar_segment_size_deg,p.seq(ss,1));
        coords_tmp = coords_tmp(abs(coords_tmp(:,1)) <= p.bar_segment_size_deg/2 & abs(coords_tmp(:,2)) <= p.seq(ss,1)/2,:);
        
        bar_dot_coords{bb} = coords_tmp(1:ndots,:);
        bar_dot_age{bb} = floor(rand(ndots,1)*p.dot_life);
        
        
        clear coords_tmp;
    end
    clear tmpndots;

    for step_num = 1:p.n_steps
        
        % ------ init step -------
        % choose directions for each dot based on seg #, which is
        % target on this segment
        %
        
        p.step_start(trial_counter) = GetSecs;
        if step_num == 1
            p.sweep_start(ss) = p.step_start(trial_counter);
        end
        
        % start w/ target segment
        for bb = 1:length(p.target_segments)
            bar_dot_dir{p.target_segments(bb)} = nan(ndots,1);
            tmp_ncoh = ceil(ndots*p.target_coherence(trial_counter));
            bar_dot_dir{p.target_segments(bb)}(1:tmp_ncoh) = 90 - 90 * p.segment_dirs(trial_counter,p.target_segments(bb));
            bar_dot_dir{p.target_segments(bb)}((tmp_ncoh+1):end) = linspace(360/(ndots-tmp_ncoh),360,ndots-tmp_ncoh);
            clear tmp_ncoh;
        end
        
        for bb = 1:length(p.sample_segments)
            bar_dot_dir{p.sample_segments(bb)} = nan(ndots,1);
            tmp_ncoh = ceil(ndots*p.sample_coherence(trial_counter));
            bar_dot_dir{p.sample_segments(bb)}(1:tmp_ncoh) = 90 - 90 * p.segment_dirs(trial_counter,p.sample_segments(bb));
            bar_dot_dir{p.sample_segments(bb)}((tmp_ncoh+1):end) = linspace(360/(ndots-tmp_ncoh),360,ndots-tmp_ncoh);
            clear tmp_ncoh;
        end
        
        
        % segment position
        if p.seq(ss,2) == 1 || p.seq(ss,2) == 2 % x is p.bar_segment_centers_deg, y is this_bar_pos(step_num)
            this_seg_pos = [p.bar_segment_centers_deg.' this_bar_pos(step_num)*ones(p.n_segments,1)];
            p.bar_pos(trial_counter,:) = [0 this_bar_pos(step_num)];
            % and update bar_pos variable
        else        %  y is p.bar_segment_centers_deg, x is this_bar_pos(step_num)
            this_seg_pos = [this_bar_pos(step_num)*ones(p.n_segments,1) (p.bar_segment_centers_deg(end:-1:1)).'];
            p.bar_pos(trial_counter,:) = [this_bar_pos(step_num) 0];
        end
        
        
        
        
        resp_yet = 0;
        % note: this won't be general for non-unitary #'s of sample
        % segments...
        %corr_resp = find(p.segment_dirs(trial_counter,p.target_segments)==p.segment_dirs(trial_counter,p.sample_segments)); % on this segment, which response button is correct?
        %fprintf('Trial %i:\tsample: %i\ttargets: %i, %i\tcorrect response: %i\n',trial_counter,p.segment_dirs(trial_counter,2),p.segment_dirs(trial_counter,1),p.segment_dirs(trial_counter,3),p.corr_resp(trial_counter));
        
        % loop over frames
        
        while (GetSecs - p.expt_start)  < (p.step_dur*trial_counter + p.start_wait)
            
            % ------ Draw Dots --------
            for bb = 1:p.n_segments
                % NOTE: when plotting, flip y - so that + is up; same for
                % centers of segments (hence the [1;-1])
                if p.seq(ss,2)==1 || p.seq(ss,2) == 2  % bar is horizontal
                    %Screen('DrawDots',w,[1; -1].*(p.ppd*bar_dot_coords{bb}(:,[1 2]).'),p.dot_size_deg*p.ppd, p.dot_color, p.ppd*this_seg_pos(bb,[1 2]) .* [1 -1] + p.scr_center, 1 );
                    Screen('DrawDots',w,repmat([1; -1],1,size(bar_dot_coords{bb},1)).*(p.ppd*bar_dot_coords{bb}(:,[1 2]).'),p.dot_size_deg*p.ppd, p.dot_color, p.ppd*this_seg_pos(bb,[1 2]) .* [1 -1] + p.scr_center, 1 );
                else  % bar is vertical, need to flip coords for dots and center
                    %Screen('DrawDots',w,[1; -1].*(p.ppd*bar_dot_coords{bb}(:,[2 1]).'),p.dot_size_deg*p.ppd, p.dot_color, p.ppd*this_seg_pos(bb,[1 2]) .* [1 -1] + p.scr_center, 1 );
                    Screen('DrawDots',w,repmat([1; -1],1,size(bar_dot_coords{bb},1)).*(p.ppd*bar_dot_coords{bb}(:,[2 1]).'),p.dot_size_deg*p.ppd, p.dot_color, p.ppd*this_seg_pos(bb,[1 2]) .* [1 -1] + p.scr_center, 1 );
                end
                
                
            end
            
            
            
            % ----- Check For Response
            % maybe wait ~500 ms to start checking? just in case continued
            % button press from previous trial?
            
            [resp, time_stamp] = checkForResp(p.resp_keys, p.esc_key);
            if resp ~= 0
                if resp == -1
                    sca;
                    fprintf('ESC pressesd during bar sweep %i, trial number %i\n',ss,trial_counter);
                    ShowCursor;
                    if p.do_et == 1
                    Eyelink('StopRecording');
                    Eyelink('ShutDown');
                    end
                    return;
                elseif resp_yet == 0
                    
                    p.rt(trial_counter) = time_stamp-p.step_start(trial_counter);
                    
                    this_resp = find(p.resp_keys==resp);
                    if (p.seq(ss,2)==3)||(p.seq(ss,2)==4)
                        this_resp = mod(this_resp,2)+1;
                    end
                    
                    p.resp(trial_counter) = this_resp;%find(p.resp_keys==resp);
                    p.correct(trial_counter) = this_resp==p.corr_resp(trial_counter);%p.resp(trial_counter)==p.corr_resp(trial_counter);
                    
                    resp_yet = 1;
                    clear this_resp
                end
            end

            
            
            
            % ----- if response, update fixation point; staircase --------
            
            % draw fixation aperture
            Screen('FillOval',w,p.bg_color,fix_aperture_rect);
            
            
            
            % draw fixation cross [[[[ for now, no updating...]]]
            if isnan(p.correct(trial_counter))
                this_fix_color = p.fix_color;
            else
                this_fix_color = p.fb_colors(~p.correct(trial_counter)+1,:);
            end
            Screen('DrawLines',w, [-1 1 0 0; 0 0 -1 1]*p.ppd*p.fix_size,p.fix_width,this_fix_color,p.scr_center);
            Screen('Flip',w);
            
            
            
            % ---- update dot positions ------
            
            for bb = 1:p.n_segments
                
                % find dead dots
                dead_dots = bar_dot_age{bb}==p.dot_life;
                
                
                % give them new starting coords (just rand(1) *
                % width/height)
                
                %bar_dot_coords{bb}(dead_dots,:) = rand(sum(dead_dots),2).*[p.bar_segment_size_deg p.seq(ss,1)] - [p.bar_segment_size_deg p.seq(ss,1)]/2;
                bar_dot_coords{bb}(dead_dots,:) = rand(sum(dead_dots),2).*repmat([p.bar_segment_size_deg p.seq(ss,1)],sum(dead_dots),1) - repmat([p.bar_segment_size_deg p.seq(ss,1)]/2,sum(dead_dots),1);
                
                % reset their age
                bar_dot_age{bb}(dead_dots) = 0;
                
                
                
                % now update all the alive dots
                dxy = [cosd(bar_dot_dir{bb}) sind(bar_dot_dir{bb})] * p.dot_speed / p.refresh_rate;
                
                
                bar_dot_coords{bb}(~dead_dots,:) = bar_dot_coords{bb}(~dead_dots,:)+dxy(~dead_dots,:);
                
                bar_dot_age{bb}(~dead_dots) = bar_dot_age{bb}(~dead_dots)+1;
                
                % find dots outside the aperture & move them
                
                tmp_xout = abs(bar_dot_coords{bb}(:,1))>(p.bar_segment_size_deg/2);
                bar_dot_coords{bb}(tmp_xout,1) = bar_dot_coords{bb}(tmp_xout,1)*-1;
                
                tmp_yout = abs(bar_dot_coords{bb}(:,2))>(p.seq(ss,1)/2);
                bar_dot_coords{bb}(tmp_yout,2) = bar_dot_coords{bb}(tmp_yout,2)*-1;
                clear dxy dead_dots tmp_xout tmp_yout;
            end
            
           
        end
        
        
        % ------- update coherence w/ staircase --------
        
        % if incorrect_counter == n_incorrect on this trial, increase coherence on next trial,
        % reset incorrect_counter
        
        % figure out staircase stuff for next trial;
        if p.correct(trial_counter) == 1
            up_cnt = up_cnt+1;
        else
            down_cnt = down_cnt+1;
        end
        
        next_coh = p.target_coherence(trial_counter); % unless it changes
        
        if up_cnt==p.up_thresh
            up_cnt = 0; down_cnt = 0;
            % update coherence  (increase difficulty, so smaller coherence
            if p.target_coherence(trial_counter)>p.coh_step % can only step if possible
                next_coh = max(p.target_coherence(trial_counter)-p.coh_step,p.coh_step); % minimum coherence allowed is the step val, don't use 0%
            end
        end
        
        if down_cnt==p.down_thresh
            down_cnt = 0; up_cnt = 0;
            % update sf step
            if p.target_coherence(trial_counter)<1 % make sure we don't jump outside our range
                next_coh = min(p.target_coherence(trial_counter)+p.coh_step,1); % easier, so bigger coherence
            end
        end

        
        
        trial_counter = trial_counter+1;
        
        if trial_counter <= (p.n_steps * size(p.seq,1))
            p.target_coherence(trial_counter) = next_coh;
        else % if this is the last trial, save this out
            tmpfn = sprintf('./data/%s_r%02.f_lastCoh.mat',p.subj,p.run);
            init_coh = next_coh;
            save(tmpfn,'init_coh');
            clear init_coh tmpfn;
        end
        
        clear this_fix_color; % and other vars
        
    end
    p.sweep_end(ss)=GetSecs;
    
end


Screen('DrawLines',w, [-1 1 0 0; 0 0 -1 1]*p.ppd*p.fix_size,p.fix_width,p.fix_color,p.scr_center);
Screen('Flip',w);

% END OF EXPERIMENT
if p.do_et == 1
    Eyelink('Message','xDAT %i',111);
end

%end_bars = GetSecs;
% ------ end wait time (check for esc) -------------

resp = 0;
fprintf('waiting for end\n');
while (GetSecs-p.sweep_end(ss)) < p.end_wait
    [resp, ~] = checkForResp(p.start_key, p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressed during post-task wait time\n');
        ShowCursor;
        if p.do_et==1
        Eyelink('StopRecording');
        Eyelink('ShutDown');
        end
        return;
    end
end
p.expt_end = GetSecs;

clear resp;


save(p.filename,'p');



% ------ end of run - present feedback -------------
Screen('Preference','TextRenderer',1);
Screen('DrawLines',w, [-1 1 0 0; 0 0 -1 1]*p.ppd*p.fix_size,p.fix_width,p.fix_color,p.scr_center);
txt = sprintf('Acc: %.02f%%, %i missed responses',100*nansum(p.correct)/sum(~isnan(p.correct)),sum(isnan(p.correct)));
Screen('TextSize', w, 30);
DrawFormattedText(w,txt,'center',p.scr_center(2)-3*p.ppd,p.fix_color);
%normBoundsRect = Screen('TextBounds', w, txt);
%txtloc = [p.xc - normBoundsRect(3)/2, 2 * p.ppd + p.yc - normBoundsRect(4)/2];
%Screen('DrawText', w, txt,txtloc(1),txtloc(2), 255 );
Screen('Flip',w);

resp = 0;
fprintf('waiting for final space\n');
while resp == 0
    [resp, ~] = checkForResp([p.start_key KbName('space')],  p.esc_key);
    if resp == -1
        sca;
        fprintf('ESC pressesd during feedback time\n');
        ShowCursor;
        return;
    end
end
clear resp;


if p.do_et == 1
    Eyelink('StopRecording');
    Eyelink('ReceiveFile',[p.eyedatafile '.edf'],[p.eyedatafile '.edf']);
    
    p.eyedatafile_renamed = [p.filename(1:(end-3)) 'edf'];
    movefile([p.eyedatafile '.edf'],p.eyedatafile_renamed);
    
    Eyelink('ShutDown');
end


Screen('CloseAll');

return;