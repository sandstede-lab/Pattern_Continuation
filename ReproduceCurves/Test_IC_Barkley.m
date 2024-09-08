% test model

clear all;
close all;

path(path,'..');

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

modelpar = struct;

modelpar.model = 'Barkley';
modelpar.xdim = 2;
modelpar.threshold = 0.7;


featpar = struct;
featpar.feature = 'retract';
modelpar.sets = 'pos';
featpar.N = 1;
featpar.alpha = 1;

reproduce = 'meander2';

contpar = struct;
contpar.step_size = 0.02;
contpar.max_step_arc_length = 0.02;
contpar.min_step_arc_length = contpar.max_step_arc_length;
contpar.L = .5;
contpar.max_sim = 1;

line = struct; line.start = [0.2,0.02]; line.end = [0.7,0.02];


switch reproduce
    case 'retract'
        % Find retract
        featpar.feature = 'retract';
        start.point = [0.375,0.015];
        start.normal = [0.01,0.01];
        contpar.step_size = 0.02;
        contpar.max_step_arc_length = 0.02;
        contpar.min_step_arc_length = 0.02;
        line = struct; line.start = [0.2,0.02]; line.end = [0.4,0.02];
        contpar.stopping = 0.9;
        modelpar.sets = 'pos';
        featpar.avoid_steady = 0;
    case 'meander1'
        start.point = [0.40,0.02];%+0.02;
        start.normal = [0.01,0.01];
        featpar.feature = 'meander';
        modelpar.sets = 'tip_points';
        line = struct; line.start = [0.3,0.02]; line.end = [0.6,0.02];
        contpar.stopping = 0.5;
        featpar.avoid_steady = 1;
    case 'meander2'
        start.point = [0.66,0.01];%+0.02;
        start.normal = [0.01,0.01];
        featpar.feature = 'meander';
        modelpar.sets = 'tip_points';
        line = struct; line.start = [0.6,0.02]; line.end = [0.8,0.02];
        contpar.stopping = 0.5;
        featpar.avoid_steady = 1;
    case 'drift'
        start.point = [0.45,0.005];
        start.normal = [0.01,0.01];
        % contpar.max_step_arc_length = 0.005*2;%*1.5;
        % contpar.min_step_arc_length = 0.005*2;%*1.5;
        % contpar.step_size = 0.005*2;%*1.5;
        contpar.stopping = 0.5;
        featpar.feature = 'drift';
        modelpar.sets = 'tip_points';
        line = struct; line.start = [0.42,0.02]; line.end = [0.6,0.02];
        featpar.avoid_steady = 1;




end

%line1 = struct; line1.start = [0.2,0.015]; line1.end = [0.4,0.015];

angle = struct;
angle.start = 0; angle.end = 9/10*pi;

cd ..
[start, dists_c1, dists_c2] = find_initial_condition(line, angle, contpar.step_size, modelpar, featpar,feature_handle, obj_handle);
cd ReproduceCurves


% start.point = [0.4959    0.0200];
% start.normal = [0.9397 0.3420];

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

cd ..
[p_history,counts,L_history,p_history_all,metric_history_all] = continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start )
cd ReproduceCurves