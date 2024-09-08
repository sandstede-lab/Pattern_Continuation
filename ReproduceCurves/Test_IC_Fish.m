% Test finding initial points

% test model

% clear all;
% close all;

rng(5);

path(path,'..');

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

modelpar = struct;
modelpar.model = 'Bullara';
modelpar.init = 'random';
modelpar.sets = 'homogeneous';
modelpar.xdim = 2;
modelpar.threshold = 0.01;

featpar = struct;


contpar = struct;
contpar.step_size = 0.1;
contpar.max_step_arc_length = 0.1;
contpar.min_step_arc_length = contpar.max_step_arc_length;%/2;
contpar.L = 20;
contpar.max_sim = 5;

% reproduce = 'steady';
% reproduce = 'neg-spot';
% reproduce = 'spot'
switch reproduce
    case 'steady'
    

        featpar.feature = 'steady';
        contpar.max_sim = 3;
        contpar.stopping = .01;
        featpar.avoid_steady = 0;
        featpar.alpha = .8;
        featpar.N = 1%contpar.max_sim;
        contpar.step_size = 0.1%25;
        contpar.max_step_arc_length = contpar.step_size;
        contpar.min_step_arc_length = contpar.max_step_arc_length;%/2;
        contpar.L = 30;
        line = struct; line.start = [0,3]; line.end = [1.8,3];
        
        line = struct; line.start = [0,20]; line.end = [1.8,20];
        
        angle = struct; angle.start = -pi; angle.end = 0*pi;
    case 'spot'
        featpar.feature = 'num'%roundness-bag';
        contpar.max_sim = 15%20;
        contpar.stopping = 1;%
        % contpar.stopping = .01;
        featpar.avoid_steady = 1;
        modelpar.sets = 'pos'%-filter';
        modelpar.threshold = 0.7;
        featpar.alpha = 2;%2%3;
        
        line = struct; line.start = [1.5,20]; line.end = [2,20];
        % line = struct; line.start = [1.5,6]; line.end = [2.2,6];
        
        featpar.N = contpar.max_sim;
        contpar.step_size = 0.3%25;
        contpar.max_step_arc_length = contpar.step_size;
        contpar.min_step_arc_length = contpar.max_step_arc_length;%/2;
        % line = struct; line.start = [0,18]; line.end = [1.8,18.1];
        
        angle = struct; angle.start = -pi; angle.end = 0*pi;
    case 'neg-spot'
        featpar.feature = 'num'%roundness-bag';
        contpar.max_sim = 15;
        contpar.stopping = 1;%.01;
        % contpar.stopping = .01;
        featpar.avoid_steady = 1;
        
        modelpar.sets = 'neg';
        % modelpar.sets = 'neg-filter';
        modelpar.threshold = 0.3;
        featpar.alpha = 2;
        featpar.N = contpar.max_sim;
        line = struct; line.start = [2.5,20]; line.end = [3,20];
        % line = struct; line.start = [2.5,6]; line.end = [3,6];
        contpar.step_size = .3%0.3;
        contpar.max_step_arc_length = contpar.step_size;
        contpar.min_step_arc_length = contpar.max_step_arc_length;%/2;
        %line = struct; line.start = [0,18]; line.end = [1.8,18.1];
        
        angle = struct; angle.start = -pi; angle.end = 0*pi;


end
%num_search = 10;
cd ..
[start, dists_c1, dists_c2] = find_initial_condition(line, angle, contpar.step_size, modelpar, featpar,feature_handle, obj_handle);
% cd ReproduceCurves
% 
% 
% 
% cd ..
[p_history,counts,L_history,p_history_all,metric_history_all] = continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start )
cd ReproduceCurves