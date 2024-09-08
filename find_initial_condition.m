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

    


% featpar.N = 1;

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
while norm(p_right - p_left) >= .1*arclength
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

start.point = p_mid;

% Then search for the optimal direction


% angles = linspace(angle.start,angle.end,num_search);
angle_left = angle.start;
angle_right = angle.end;

modelpar.a = start.point(1) + arclength*cos(angle_left);
modelpar.b = start.point(2) + arclength*sin(angle_left);
featuresl = feature_handle(featpar, modelpar);
if iscell(featuresl)
    featuresl = featuresl{1};
end

modelpar.a = start.point(1) + arclength*cos(angle_right);
modelpar.b = start.point(2) + arclength*sin(angle_right);
featuresr = feature_handle(featpar, modelpar);
if iscell(featuresr)
    featuresr = featuresr{1};
end

search_2 = [];
dists_c2 = [];



while norm(angle_right - angle_left) >= pi/8

    angle_mid = 0.5*( angle_left + angle_right );
    modelpar.a = start.point(1) + arclength*cos(angle_mid);
    modelpar.b = start.point(2) + arclength*sin(angle_mid);
    featuresm = feature_handle(featpar, modelpar);
    if iscell(featuresm)
        featuresm = featuresm{1};
    end
    obj_l = obj_handle(featuresm,featuresl);
    search_2 = [search_2, start.point'+arclength*[cos(0.5*(angle_left + angle_mid)');sin(0.5*(angle_left + angle_mid))]  ];
    dists_c2 = [dists_c2, obj_l];

    obj_r = obj_handle(featuresm,featuresr);
    search_2 = [search_2, start.point'+arclength*[cos(0.5*(angle_right + angle_mid)');sin(0.5*(angle_right + angle_mid))]  ];
    dists_c2 = [dists_c2, obj_r];

    
    if obj_l > obj_r
        angle_right = angle_mid;
        featuresr = featuresm;
    else
        angle_left = angle_mid;
        featuresl = featuresm;
    end



end

angle_mid = 0.5*( angle_left + angle_right );
%idx = find(dists_c2 == max(dists_c2)); %idx = idx( round(length(idx)/2));
start.normal = [cos(angle_mid),sin(angle_mid)];
    

figure(19);
scatter( search_1(1,:), search_1(2,:), [], dists_c1 , 'filled' )
hold on,
scatter( search_2(1,:), search_2(2,:), [], dists_c2 , 'filled'  )
        

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
        xlim([0.035,0.065])
        ylim([0.005,0.06])
        ylabel('f')
        xlabel('k')

    case 'Schnakenberg'
        xlim([4.5,7.5]);
        ylim([4.5,6.5]);
        xlabel('a')
        ylabel('b')
    case 'Bullara'
        xlim([0,4])
        ylim([0,20])
        xlabel('l_x')
        ylabel('h')
end
set(gca,'fontsize',24);