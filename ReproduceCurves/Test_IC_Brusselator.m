% Test finding initial points

% test model

clear all;
close all;

path(path,'..');

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

modelpar = struct;
modelpar.model = 'Brusselator';
modelpar.init = 'random';
modelpar.sets = 'pos';
modelpar.xdim = 2;
modelpar.threshold = 0.7;

featpar = struct;
featpar.feature = 'num';
featpar.alpha = 1;


contpar = struct;
contpar.step_size = 0.5;
contpar.max_step_arc_length = 0.5;
contpar.min_step_arc_length = contpar.max_step_arc_length/2;
contpar.L = 8;
contpar.max_sim = 10;

reproduce = 'neg';

switch reproduce
    case 'pos'
        
        modelpar.sets = 'pos';
        modelpar.threshold = 0.7;
        contpar.stopping = 10;
        featpar.avoid_steady = 1;
    case 'neg'

        start.point = [4.2,6.5];
        start.normal = [0,1];
        modelpar.sets = 'neg';
        modelpar.threshold = 0.3;
        contpar.stopping = 10;
        featpar.avoid_steady = 1;

    case 'steady'
        start.point = [4.5,6.5];
        start.normal = [0,1];
        featpar.feature = 'steady';
        contpar.max_sim = 1;
        contpar.stopping = .1;
        featpar.avoid_steady = 0;

end

line = struct; line.start = [3,7]; line.end = [6,7];
featpar.N = 3;
angle = struct; angle.start = 0; angle.end = pi*9/10;
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