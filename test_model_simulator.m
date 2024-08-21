% test model

clear all;
close all;

modelpar = struct;

modelpar.model = 'Brusselator';
modelpar.a = 3.5;
modelpar.b = 6;
modelpar.sets = 'pos';
modelpar.xdim = 2;
modelpar.threshold = 0.7;

% pattern = model(modelpar);
% 
% figure;
% scatter(pattern(:,1), pattern(:,2));

featpar = struct;
featpar.feature = 'num';
featpar.N = 3;
featpar.alpha = 1;

% features1 = feature_evaluation(featpar, modelpar);

%figure;histogram(features1);

% modelpar.a = 3.5%4.0;
% modelpar.b = 6;
% features2 = feature_evaluation(featpar, modelpar);
% 
% objective_evaluation(features1, features2)


start.point = [3.5,6];
start.normal = [0,1];
modelpar.sets = 'pos';
modelpar.threshold = 0.7;

contpar = struct;
contpar.step_size = 0.5;
contpar.max_step_arc_length = 0.5;
contpar.min_step_arc_length = contpar.max_step_arc_length/2;
contpar.L = 8;

contpar.max_sim = 5;

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

[p_history,counts,L_history,p_history_all,metric_history_all] = continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start )
