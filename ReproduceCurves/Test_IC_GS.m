% Test finding initial points

% test model

rng(5)
clear all;
close all;

path(path,'..');

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

modelpar = struct;
modelpar.model = 'GS';
modelpar.init = 'random';
modelpar.sets = 'pos';
modelpar.xdim = 2;
modelpar.threshold = 0.7;

featpar = struct;
featpar.feature = 'num';
featpar.alpha = 0.01;


contpar = struct;
contpar.step_size = 0.001/2;
contpar.max_step_arc_length = 0.001/2;
contpar.min_step_arc_length = contpar.max_step_arc_length;
contpar.L = 0.1;
contpar.max_sim = 10;

reproduce = 'pos';

switch reproduce
    case 'pos'
        featpar.feature = 'roundness-bag';
        modelpar.sets = 'pos';
        modelpar.threshold = 0.7;
        contpar.stopping = .01;
        %init1 = [0.053,0.024]; init2 = [0.053,0.025];
        featpar.avoid_steady = 0;
        % contpar.stopping = 0.1;
        line = struct; line.start = [0.052,0.024]; line.end = [0.052,0.023];
        % line.start = [0.06,0.05]; line.end = [0.065,0.05];
        
        % line = struct; line.start = [0.053,0.026]; line.end = [0.053,0.024];
        % line = struct; line.start = [0.051,0.015]; line.end = [0.051,0.022];

        %line = struct; line.start = [0.049,0.01]; line.end = [0.049,0.02];


        
    case 'neg'

        featpar.feature = 'roundness-bag';
        modelpar.sets = 'neg';
        modelpar.threshold = 0.3;
        contpar.stopping = .01%*10;
        %init1 = [0.053,0.024]; init2 = [0.053,0.025];
        featpar.avoid_steady = 0;
        % contpar.stopping = 0.1;
        line = struct; line.start = [0.053,0.021]; line.end = [0.053,0.024];

    case 'steady-spots'
        
        modelpar.sets = 'pos';
        featpar.feature = 'num-spots';

        contpar.max_sim = 1;
        contpar.stopping = .1;
        featpar.avoid_steady = 0;
        
        line = struct; line.start = [0.053,0.024]; line.end = [0.053,0.023];

        contpar.step_size = 0.002;
        contpar.max_step_arc_length = contpar.step_size;
        contpar.min_step_arc_length = contpar.step_size;

    case 'steady'
        
        modelpar.sets = 'homogeneous';
        
        featpar.feature = 'steady';

        

        contpar.max_sim = 1;
        contpar.stopping = .1;
        featpar.avoid_steady = 0;
        line = struct; line.start = [0.045,0.02]; line.end = [0.045,0.03];
        line = struct; line.start = [0.038,0.02]; line.end = [0.038,0.015];

        line = struct; line.start = [0.04,0.02]; line.end = [0.04,0.015];

        contpar.step_size = 0.002;
        contpar.max_step_arc_length = contpar.step_size;
        contpar.min_step_arc_length = contpar.step_size;
        
        % load('steady_fal.mat')
        % start = struct;
        % start.point = p_history(:,end-3)'-0.001*[1,.5];
        % start.normal = p_history(:,end-2)' - p_history(:,end-3)';%+0.001*[0,100];
        % line = struct; line.start = [0.05,0.02]; line.end = [0.05,0.03];


end

%line = struct; line.start = [0.05,0.005]; line.end = [0.05,0.02];
featpar.N = 1;
angle = struct; angle.start = -pi/2; angle.end = pi/2;
%num_search = 10;
cd ..
% featpar.feature = 'steady';
[start, dists_c1, dists_c2] = find_initial_condition(line, angle, contpar.step_size, modelpar, featpar,feature_handle, obj_handle);
cd ReproduceCurves
% % 
% % 
% % 
% load GS_start_steady
% featpar.feature = 'roundness-bag';
[p_history,counts,L_history,p_history_all,metric_history_all] = continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start )
%cd ReproduceCurves