% test model

clear all;
close all;

modelpar = struct;

modelpar.model = 'Barkley';


modelpar.a = 0.37;
modelpar.b = 0.01;
modelpar.sets = 'tipping_points';
modelpar.xdim = 2;
modelpar.threshold = 0.7;

% pattern = model(modelpar);
% 
% figure;
% scatter(pattern(:,1), pattern(:,2));
% 
featpar = struct;
featpar.feature = 'retract';
modelpar.sets = 'pos';
featpar.N = 1;
featpar.alpha = 1;

% features1 = feature_evaluation(featpar, modelpar);

%figure;histogram(features1);

modelpar.a = 0.38%4.0;
modelpar.b = 0.02;
% features2 = feature_evaluation(featpar, modelpar);
% 
% objective_evaluation(features1, features2)

% Find retract
featpar.feature = 'retract';
start.point = [0.375,0.015];
start.normal = [0.01,0.01];


% Find meander

start.point = [0.66,0.01]%+0.02;
featpar.feature = 'meander';
modelpar.sets = 'tipping_points';

contpar = struct;
contpar.step_size = 0.01;
contpar.max_step_arc_length = 0.01;
contpar.min_step_arc_length = contpar.max_step_arc_length/2;
contpar.L = .5;


% Find drift
start.point = [0.45,0.005];
contpar.max_step_arc_length = 0.005*1.5;
contpar.min_step_arc_length = 0.005*1.5;
contpar.step_size = 0.005*1.5;
featpar.feature = 'drift';

% Switch to Rossler
featpar.feature = 'area';
modelpar.sets = 'pos';
modelpar.threshold = 0.9;
contpar.L = 1.5;
start.point = [3.9,0.17];
start.normal = [-0.1,0];
contpar.max_step_arc_length = 0.03;
contpar.min_step_arc_length = 0.03;
contpar.step_size = 0.03;
modelpar.model = 'Rossler';

% Switch to BE
contpar.L = 0.5;
modelpar.model = 'Bar-Eiswirth';

start.point = [0.07,0.03];
start.normal = [0.0,0.01];
contpar.max_step_arc_length = 0.005;
contpar.min_step_arc_length = 0.005;
contpar.step_size = 0.005;
featpar.feature = 'turbulence';
modelpar.sets = 'tipping_points';



contpar.max_sim = 1;

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start )
