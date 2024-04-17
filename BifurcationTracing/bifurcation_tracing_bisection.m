function [p_history,counts,L_history,p_history_all,metric_history_all] = bifurcation_tracing_bisection(init1, init2, L,max_step_arc_length, min_step_arc_length, dist, max_bisect, model, sets,steady,interp)

% Bifurcation tracing `algorithm' 
% N: number of steps
% dist: distance metric, can be:
% 1. 'num': number of alpha shapes, 1 number per simulation
% 2. 'perimeter': perimeter of alpha shapes, 1 number per cluster in each
% simulation
% 3. 'area': area of alpha shapes, 1 number per cluster in each
% simulation
% 4. 'roundness': roundness of alpha shapes, 1 number per cluster in each
% simulation
    

rng(5);

if nargin < 5
    interp = false;
end
if nargin < 6
    adapt = false;

    
end
if nargin < 7
    model='Brusselator';
end
if nargin < 8
    sets = 'pos';
end


if nargin < 11
    interp = 0;
end

if nargin < 12
    Uinit = [];
end



a(1) = init1(1);
b(1) = init1(2);

a(2) = init2(1);
b(2) = init2(2);

alpha = [];
max_count = 1;


% 2 initial points for transition curve

switch model
    case 'Brusselator'
        
        
        

        max_count = 10;
        
        alpha = 1;
%         alpha = 1;

        
    case {'SH','SH_spots_init','SH_stripes_init'}
        
        
        
        
        alpha = 0.3+0.05;
        max_count = 5;
        
    case 'GS'
        
        alpha = 0.01;
        max_count = 5;
        
    
    case 'Schnakenberg'
        alpha = 0.1;
        max_count = 10;
    case 'SH_1D'
        alpha = [];
        max_count = 1;
   
end

if steady  || strcmp(model,'SH_1D') || strcmp(model,'SH_spots_init') || strcmp(model,'SH_stripes_init')
    max_count = 1;
end
% For each combination of parameter, do 5 simulations




p_current = [a(2),b(2)];

p_prior = [a(1),b(1)];

p_history(:,1) = p_prior ;
p_history(:,2) = p_current ;

direction = p_current - p_prior;
direction = direction/norm(direction);
direction_norm = [-direction(2), direction(1)];
% Kept alpha as constant 1.2


% figure(1);
% 
% plot(a(1:2), b(1:2), 'ko-')
t_max = 100;

L_current = 0;
i = 0;



p_prior = p_current;
    


step_size = min_step_arc_length;

% 
flag = 0;

flag0 = flag;
p_history_all = [];
metric_history_all = [];
while L_current < L
    i = i+1;
    fprintf(['step ',num2str(i), '\n']);
    p_prior = p_current;

    p_pred = p_prior + step_size*direction;
    
    
    % Tentative correction of prediction along normal direction
    p_c1 = p_pred + 0.5*step_size*direction_norm;
    p_c2 = p_pred - 0.5*step_size*direction_norm;
    
    
    
    dn = .5*step_size*direction_norm;    
            
    
    
    
    [p_current,counts(i),p_history_new, metric_history_new] = bisect_significant_interval(model, p_c1, p_c2, dn, max_bisect, alpha,max_count, dist,t_max,sets,steady, interp);
    p_history_all = [p_history_all, p_history_new];
    metric_history_all = [metric_history_all, metric_history_new];
    p_history(:,i+2) = p_current';
    L_history(i) = step_size;
    L_current = L_current + step_size;%sqrt(sum( (p_current - p_prior).^2 ));
    
    if counts(i) <= 2
        step_size = max( [ step_size/sqrt(sqrt(2)),min_step_arc_length ]);
        fprintf(['step_size decreased to',num2str(step_size),'\n']);
    end
    if counts(i) >= max_count
        step_size = min([sqrt(sqrt(2))*step_size, max_step_arc_length]) ;
        fprintf(['step_size increased to',num2str(step_size),'\n']);
    end
%     
    
  
    
    % Update direction vectors
    direction = p_current - p_prior;
    direction = direction/norm(direction);
    direction_norm = [-direction(2), direction(1)];
    
    
    if i>5 && (  (max( metric_history_new ) <0.5) || (min( metric_history_new >-0.5) )) && strcmp(dist,'line_fit')
        L_current = 2*L;
    end
    if i>5 && ( ( max( metric_history_new ) <0.5 )||( min( metric_history_new) >0.5 )) && strcmp(dist,'ring_fit')
        L_current = 2*L;
    end
    if i>5 && ( ( max( metric_history_new ) <0.5 )||( min( metric_history_new) >0.5 )) && ( strcmp(dist,'steady') || strcmp(dist,'turbulence')  || strcmp(dist,'turbulence2') )
        L_current = 2*L;
    end
    
    
    
    
    
        figure(19)
    %     hold on,

        plot(p_history(1,1:i+2), p_history(2,1:i+2), 'ko-')
        hold on,
        scatter(p_history_all(1,:), p_history_all(2,:),[],metric_history_all,'filled');
        
        
        hold off
        switch model
            case 'Brusselator'

                ylim([6,15])
                xlim([3.5,9])
                xlabel('a')
                ylabel('b')
            case {'SH','SH_stripes_init','SH_spots_init'}
                ylim([0,2])
                xlim([-0.15,0.4])
                xlabel('\mu')
                ylabel('\nu')
                
            case 'GS'
                xlim([0.05,0.065])
                ylim([0.023,0.038])
                ylabel('f')
                xlabel('k')
                
            case 'SH_1D'
                xlim([0,1.2]);
                ylim([1.7,3]);
                xlabel('\mu')
                ylabel('\nu')
            case 'Schnakenberg'
                xlim([4.5,13]);
                ylim([4.5,7.5]);
                xlabel('a')
                ylabel('b')
                
                
            case 'Barkley'

                ylim([0.01,0.13])
                xlim([0.25,0.9])
                xlabel('a')
                ylabel('b')
            case 'Bar-Eiswirth'
                xlim([0.0,0.15])
                ylim([0.04,0.26])
                xlabel('a')
                ylabel('\epsilon')
            case 'Rossler'
                xlim([2.5,4])
                ylim([0.16,0.26])
                xlabel('c')
                ylabel('a')

        end
        end
        title(['i=',num2str(i)]);
    
    pause(0.1)
%     
    
    
    
   
    
    
    
end
