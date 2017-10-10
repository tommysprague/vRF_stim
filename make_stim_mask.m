% function make_stim_mask.m
%
% creates stimulus mask using p struct "p", saves to file "fn_out", renders
% at resolution "res" (default 270x270)
%
% fn_out should just be a prefix, no .mat - will make an _images and
% _params file for vista
%
% NOTE: p.start_wait, end_wait should be integer # of TRs... (as should
% step_dur)
%
% TCS, 6/12/2017


function make_stim_mask(p,fn_out,TR,res)


if nargin < 4
    res = 270;
end

% if TR in ms...
if TR > 100
    TR = TR/1000;
end

dir_all = nan(size(p.bar_pos,1),1);
size_all = nan(size(p.bar_pos,1),1);
idx = 1;

for bb = 1:size(p.seq,1)
    size_all(idx:(idx+p.n_steps-1)) = p.seq(bb,1)*ones(p.n_steps,1);
    dir_all(idx:(idx+p.n_steps-1)) = p.seq(bb,2)*ones(p.n_steps,1);
    idx = idx+p.n_steps;
end

gg = linspace(min(min(p.bar_pos - size_all/2,[],2),[],1),max(max(p.bar_pos + size_all/2,[],2),[],1),res);

TRs_per_step = round(p.step_dur/TR);

n_TRs = round((p.expt_end-p.expt_start)/TR);
fprintf('In total, %i TRs\n',n_TRs);


% ------- PARAMS FILE ------
% need a _params file, which has stimulus.seq and stimulus.seqtiming
% stimulus.seq is 1:n_TRs

stimulus.seq = 1:n_TRs;

% stimulus.seqtiming is TR*(stimulus.seq-1)
stimulus.seqtiming = TR*(stimulus.seq-1);

fn_seq = sprintf('%s_params.mat',fn_out);
save(fn_seq,'stimulus');


% ------ IMAGES FILE ------

images = zeros(res,res,n_TRs);

start_TR = 1+round(p.start_wait/TR); % figure out when to start (1+TR/p.start_wait)


% keep track of where we are in iamges
this_TR = start_TR;

for tt = 1:size(p.bar_pos,1)
    
    proto_img = gg>=(sum(p.bar_pos(tt,:)) - size_all(tt)/2) & gg<=(sum(p.bar_pos(tt,:)) + size_all(tt)/2);
    
    fprintf('Step %i, %i pixels stimulated\n',tt,sum(proto_img));
    
    if dir_all(tt) < 2.5 % up/down - horizontal
        images(:,:,this_TR:(this_TR+TRs_per_step-1)) = repmat(proto_img.',1,res,TRs_per_step);
    else
        images(:,:,this_TR:(this_TR+TRs_per_step-1)) = repmat(proto_img,res,1, TRs_per_step);
    end
    
    this_TR = this_TR+TRs_per_step;
    
end

images = flipud(images); % because vista wants i,j coords


fn_imgs = sprintf('%s_images.mat',fn_out);
save(fn_imgs,'images');



return
