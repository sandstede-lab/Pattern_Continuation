function [start, dists_c1, dists_c2] = find_initial_condition(line, angle, arclength, modelpar, featpar,feature_handle, obj_handle)

% The function searches for an initial point and direction for
% continuation. 
%
% Required input:
%
% line: a struct with two attributes line.start and line.end specifying the start and end points in
% the 2D parameter space. 
%
% angle: a struct with two attributes angle.start and angle.end specifying
% the start and end of angle for the direction search. 
% 
% arclength: resolution of the search, should be consistent with
% continuation

path(path,'../Models/Reaction_Diffusion');
path(path,'../Models/Spiral_Wave');

    


featpar.N = 1;

p_left = line.start;
p_right = line.end;

modelpar.a = p_left(1);
modelpar.b = p_left(2);
featuresl = feature_handle(featpar, modelpar);
if iscell(featuresl)
    featuresl = featuresl{1};
end

modelpar.a = p_right(1);
modelpar.b = p_right(2);
featuresr = feature_handle(featpar, modelpar);
if iscell(featuresr)
    featuresr = featuresr{1};
end

search_1 = [];
dists_c1 = [];

% First do a bisection to find the starting point
while norm(p_right - p_left) > arclength
    p_mid = 0.5*( p_left + p_right );
    modelpar.a = p_mid(1);
    modelpar.b = p_mid(2);
    featuresm = feature_handle(featpar, modelpar);
    if iscell(featuresm)
        featuresm = featuresm{1};
    end
    obj_l = obj_handle(featuresm,featuresl);
    search_1 = [search_1, 0.5*(p_left + p_mid)'  ];
    dists_c1 = [dists_c1, obj_l/norm(p_left-p_mid)];

    obj_r = obj_handle(featuresm,featuresr);
    search_1 = [search_1, 0.5*(p_right + p_mid)'  ];
    dists_c1 = [dists_c1, obj_r/norm(p_right-p_mid)];

    if obj_l > obj_r
        p_right = p_mid;
        featuresr = featuresm;
    else
        p_left = p_mid;
        featuresl = featuresm;
    end



end



% Then search for the optimal direction

num_search = 10;
angles = linspace(angle.start,angle.end,num_search);

start = struct;
start.point = ( p_left + p_right )/2;

featuresl = cell(1, num_search);
featuresr = cell(1, num_search);

% Record result of dists on search 1
for counter = 0:num_search-1
    next = start.point + arclength*[cos(angles(counter+1)),sin(angles(counter+1))  ];
    normal = [sin(angles(counter+1)), -cos(angles(counter+1))];

    modelpar.a = next(1) + .5*arclength*normal(1);
    modelpar.b = next(2) + .5*arclength*normal(2);
    featuresl{counter+1} = feature_handle(featpar, modelpar);
    if iscell(featuresl{counter+1})
        featuresl{counter+1} = featuresl{counter+1}{1};
    end


    modelpar.a = next(1) - .5*arclength*normal(1);
    modelpar.b = next(2) - .5*arclength*normal(2);
    featuresr{counter+1} = feature_handle(featpar, modelpar);
    if iscell(featuresl{counter+1})
        featuresr{counter+1} = featuresr{counter+1}{1};
    end
    dists_c2(counter+1) = obj_handle(featuresl{counter+1},featuresr{counter+1})/norm(arclength*normal);




end

idx = find(dists_c2 == max(dists_c2)); %idx = idx( round(length(idx)/2));
start.normal = [cos(mean(angles(idx))), sin(mean(angles(idx)))];

    

figure(19);
scatter( search_1(1,:), search_1(2,:), [], dists_c1 , 'filled' )
hold on,
scatter( start.point(1)+cos(angles)*arclength, start.point(2)+sin(angles)*arclength, [], dists_c2 , 'filled'  )
        

switch modelpar.model
    case 'Brusselator'
        
        ylim([6,13])
        xlim([3,7])
        xlabel('a')
        ylabel('b')
    case 'SH'
        ylim([0,2])
        xlim([-0.15,0.4])
        xlabel('\mu')
        ylabel('\nu')

    case 'GS'
        xlim([0.05,0.065])
        ylim([0.023,0.038])
        ylabel('F')
        xlabel('k')


    case 'Schnakenberg'
        xlim([4.5,7.5]);
        ylim([4.5,6.5]);
        xlabel('a')
        ylabel('b')
end
set(gca,'fontsize',24);