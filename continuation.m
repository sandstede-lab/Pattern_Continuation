function [p_history,counts,L_history,p_history_all,metric_history_all] = continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start, viz )

if nargin < 7
    viz = 1;
end

start.normal = start.normal/norm(start.normal);

p_current = start.point;%+contpar.step_size*start.normal;

p_prior = start.point;

p_history(:,1) = p_prior ;
p_history(:,2) = p_current ;

direction = start.normal;
direction = direction/norm(direction);
direction_norm = [-direction(2), direction(1)];

L_current = 0;
i = 0;

max_count = contpar.max_sim;

p_prior = p_current;
    
min_step_arc_length = contpar.min_step_arc_length;
max_step_arc_length = contpar.max_step_arc_length;

step_size = contpar.step_size;

% 
flag = 0;

flag0 = flag;
p_history_all = [];%p_history;
metric_history_all = [];
L_history = [];

while L_current < contpar.L
    i = i+1;
    fprintf(['step ',num2str(i), '\n']);
    p_prior = p_current;

    p_pred = p_prior + step_size*direction;
    
    
    % Tentative correction of prediction along normal direction
    p_c1 = p_pred + 0.5*step_size*direction_norm;
    p_c2 = p_pred - 0.5*step_size*direction_norm;
    
    
    
    dn = .5*step_size*direction_norm;    
            
    
    
    
    [p_current,counts(i),p_history_new, metric_history_new, val_c1, val_c2] = bisect_interval(modelpar, featpar, feature_handle, obj_handle, p_c1, p_c2, dn, max_count);
    if max([val_c1,val_c2]) < contpar.stopping
        L_current = 2*contpar.L;
    else
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
    end
%     
    
  
    
    % Update direction vectors
    direction = p_current - p_prior;
    direction = direction/norm(direction);
    direction_norm = [-direction(2), direction(1)];
    
    'current prediction:'
    p_current

    'current direction'
    direction
    % if i>5 && (  (max( metric_history_new ) <0.5) || (min( metric_history_new >-0.5) )) && strcmp(featpar.feature,'line_fit')
    %     L_current = 2*L;
    % end
    % if i>5 && ( ( max( metric_history_new ) <0.5 )||( min( metric_history_new) >0.5 )) && strcmp(featpar.feature,'ring_fit')
    %     L_current = 2*L;
    % end
    % if i>5 && ( ( max( metric_history_new ) <0.5 )||( min( metric_history_new) >0.5 )) && ( strcmp(featpar.feature,'steady') || strcmp(featpar.feature,'turbulence')  || strcmp(featpar.feature,'turbulence2') )
    %     L_current = 2*L;
    % end

    
    
        if viz
    
    
    
            figure(19)
        %     hold on,
            if ~isempty(p_history)
                plot(p_history(1,:), p_history(2,:), 'ko-')
            end
            hold on,
            if ~isempty(metric_history_all)
                scatter(p_history_all(1,:), p_history_all(2,:),[],metric_history_all,'filled');
            end
            
            
            hold off
            switch modelpar.model
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
                    xlim([0.035,0.065])
                    ylim([0.005,0.08])
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
                case 'Bullara'
                    xlim([0,4])
                    ylim([0,20])
                    xlabel('l_x')
                    ylabel('h')
    
            end
            title(['i=',num2str(i)]);
        
            pause(0.1)
        end
    end
    
%     
    
    
    
   
    
    
    
end

% 
% function [p_current,counter, sim_p_history, sim_metric_history]=bisect_interval(modelpar, featpar, feature_handle, obj_handle, p_c1, p_c2, dn, max_count)
% 
% interp = 1;
% max_bisect = 2;
% 
% 
% counter_bisect = 0;
% flag = 0;
% counter = 0;
% feature_bp_all = [];
% feature_bm_all = []; 
% feature_b0_all = [];
% sim_p_history = [];
% sim_metric_history = [];
% 
% counter = 0;
% 
% feature_bp = [];feature_bm = [];feature_b0 = [];
% val_c2 = 0; val_c1 = 0; 
% dists_c2 = []; dists_c1 = []; 
% feature_bp_old = [];
% feature_bm_old = [];
% 
% p_c0 = p_c1 + p_c2; p_c0 = p_c0/2;
% while counter_bisect < max_bisect
% 
% 
%     flag = 0;
%     counter = 0;
%     feature_bp = [];feature_bm = [];feature_b0 = [];
%     val_c2 = 0; val_c1 = 0; 
%     dists_c2 = []; dists_c1 = []; 
% 
% 
% 
%     while flag == 0
% 
% 
%         if counter < max_count
%             if counter+1 > length(feature_bp_old) || counter_bisect == 0
%                 % switch modelpar.model
%                 %     case {'Barkley','Bar-Eiswirth','Rossler'}
%                 %         [dists_c1(counter+1),~,~,feature_bp{counter+1},feature_b0{counter+1}] = local_distance_vn_spiral_c_alpha(p_c1(1), p_c1(2), dn, 1,alpha,dist,t_max,model,[],[],sets,steady);
%                 % 
%                 %         [dists_c2(counter+1),~,~,~,feature_bm{counter+1}] = local_distance_vn_spiral_c_alpha(p_c2(1), p_c2(2), dn, 1,alpha,dist,t_max,model,feature_b0{counter+1},[],sets,steady);
%                 % 
%                 %     otherwise
%                         featpar.N = 1;
%                         modelpar.a = p_c1(1) + dn(1);
%                         modelpar.b = p_c1(2) + dn(2);
% 
% 
%                         feature_bp{counter+1} = feature_handle(featpar, modelpar);
%                         if iscell(feature_bp{counter+1})
%                             feature_bp{counter+1} = feature_bp{counter+1}{1};
%                         end
%                         modelpar.a = p_c1(1) - dn(1);%( p_c1(1) + p_c2(1) )/2;
%                         modelpar.b = p_c1(2) - dn(2);%( p_c1(2) + p_c2(2) )/2;
% 
%                         feature_b0{counter+1} = feature_handle(featpar, modelpar);
%                         if iscell(feature_b0{counter+1})
%                             feature_b0{counter+1} = feature_b0{counter+1}{1};
%                         end
% 
%                         dists_c1(counter+1) = obj_handle(feature_bp{counter+1},feature_b0{counter+1});
% 
%                         modelpar.a = p_c2(1) - dn(1);
%                         modelpar.b = p_c2(2) - dn(2);
% 
% 
%                         feature_bm{counter+1} = feature_evaluation(featpar, modelpar);
%                         if iscell(feature_bm{counter+1})
%                             feature_bm{counter+1} = feature_bm{counter+1}{1};
%                         end
% 
%                         dists_c2(counter+1) = obj_handle(feature_bm{counter+1},feature_b0{counter+1});
% 
% 
% 
%                 % end
%             else
%                 % switch modelpar.model
%                     % case {'Barkley','Bar-Eiswirth','Rossler'}
%                     %     feature_bp{counter+1} = feature_bp_old{counter+1};
%                     %     [dists_c1(counter+1),~,~,~,feature_b0{counter+1}] = local_distance_vn_spiral_c_alpha(p_c1(1), p_c1(2), dn, 1,alpha,dist,t_max,model,feature_bp{counter+1},[],sets,steady);
%                     %     feature_bm{counter+1} = feature_bm_old{counter+1};
%                     %     [dists_c2(counter+1),~,~,~,~] = local_distance_vn_spiral_c_alpha(p_c2(1), p_c2(2), dn, 1,alpha,dist,t_max,model,feature_b0{counter+1},feature_bm{counter+1},sets,steady);
%                     % 
%                     % otherwise
%                         feature_bp{counter+1} = feature_bp_old{counter+1};
%                         feature_bm{counter+1} = feature_bm_old{counter+1};
% 
%                         modelpar.a = p_c1(1) + dn(1);
%                         modelpar.b = p_c1(2) + dn(2);
% 
% 
%                         feature_bp{counter+1} = feature_handle(featpar, modelpar);
%                         if iscell(feature_bp{counter+1})
%                             feature_bp{counter+1} = feature_bp{counter+1}{1};
%                         end
% 
%                         modelpar.a = p_c1(1) - dn(1);%( p_c1(1) + p_c2(1) )/2;
%                         modelpar.b = p_c1(2) - dn(2);%( p_c1(2) + p_c2(2) )/2;
% 
% 
%                         feature_b0{counter+1} = feature_handle(featpar, modelpar);
%                         if iscell(feature_b0{counter+1})
%                             feature_b0{counter+1} = feature_b0{counter+1}{1};
%                         end
%                         dists_c1(counter+1) = obj_handle(feature_bp{counter+1},feature_b0{counter+1});
% 
%                         modelpar.a = p_c2(1) - dn(1);
%                         modelpar.b = p_c2(2) - dn(2);
% 
% 
%                         feature_bm{counter+1} = feature_evaluation(featpar, modelpar);
%                         if iscell(feature_bm{counter+1})
%                             feature_bm{counter+1} = feature_bm{counter+1}{1};
%                         end
%                         dists_c2(counter+1) = obj_handle(feature_bm{counter+1},feature_b0{counter+1});
% 
% 
%                 % end
%             end
% %         end
%             counter = counter + 1;
% %         else
% %             counter = max_count;
% 
%         end
% 
% 
% 
%         feature_bm_all = [];
%         feature_bp_all = [];
%         feature_b0_all = [];
% 
%         for ii = 1:counter
%                 switch featpar.feature
%                     case 'roundness'
%                         feature_bm_all{ii} = feature_bm{ii};
% 
%                         feature_bp_all{ii} = feature_bp{ii};
% 
%                         feature_b0_all{ii} = feature_b0{ii};
%                     otherwise
% 
%                         feature_bm_all = [feature_bm_all,feature_bm{ii}];
% 
%                         feature_bp_all = [feature_bp_all,feature_bp{ii}];
% 
%                         feature_b0_all = [feature_b0_all,feature_b0{ii}];
%                 end
% 
% 
% 
% %                 end
%         end
% 
% %         feature_bm_all{1}
% %         feature_bp_all{1}
% %         feature_b0_all{1}
% 
% %         [~,idx] = sort([mean(dists_c1(:)),mean(dists_c2(:))]);
% 
%         if counter > 1
%             flag = ttest(dists_c1(:),dists_c2(:));
%         end
%         % if counter > 1 || max_count == 1
% 
%         val_c2 = obj_handle(feature_bm_all, feature_b0_all);
%         val_c1 = obj_handle(feature_bp_all, feature_b0_all);
% 
%         % end
%         if featpar.avoid_steady == 1% avoid homogeneous unless we want that
%             if feature_bp_all == 0
%                 val_c1 = -999;
%             end
%             if feature_bm_all == 0
%                 val_c2 = -999;
%             end
%         end
% 
%         sim_p_history = [sim_p_history, p_c1'+dn', p_c1'-dn', p_c2'-dn'];
% 
%         sim_metric_history = [sim_metric_history, mean(cat(2,feature_bp_all)), mean(cat(2,feature_b0_all)), mean(cat(2,feature_bm_all))];
% 
%         if counter >= max_count
%             flag = 1;
%         end
%     end
% 
% 
%     val_c1
%     val_c2
% 
% %     if counter_bisect + 1 < max_bisect
% %         p_c2_prior = p_c2; 
% %         p_c1_prior = p_c1;
%         counter_bisect = counter_bisect + 1;
%         if counter_bisect == 1
%             val_c2_prior = val_c2;
%             val_c1_prior = val_c1;
%             p_c2_prior = p_c2; 
%             p_c1_prior = p_c1;
% 
%         end
% 
%         if val_c2 > val_c1 
%     %         if steady == 0 
%                 choice = 'left';
%                 p_current = p_c2;
%                 p_prev = p_c1;
%                 if counter_bisect  < max_bisect
% 
% 
% 
%                     if counter_bisect == 1
%                         dn = 0.5*dn;
% 
%                         p_c2 = p_c2 - dn;
%                         p_c1 = p_c2 + 2*dn;
%                         %p_c1 = p_c2 + dn;
%                     end
% 
%                 end
% 
%                 feature_bp_old = feature_b0; 
%                 feature_bm_old = feature_bm;
% 
%     %         else
%     %             
%     %         end
% 
% 
%         else
%             choice = 'right';
%             p_current = p_c1;
%             p_prev = p_c2;
% 
%             if counter_bisect  < max_bisect
% 
%                 val_c2_prior = val_c2;
%                 val_c1_prior = val_c1;
% 
%                 if counter_bisect == 1
%                     dn = 0.5*dn;
% 
%                     p_c2 = p_c1 - dn;
%                     p_c1 = p_c2 + 2*dn;
%                     %p_c1 = p_c1 + dn;
%                 end
%             end
% 
%             feature_bm_old = feature_b0; 
%             feature_bp_old = feature_bp;
%         end
% 
% %     end
% 
% %     flag = ttest(dists_c1(:),dists_c2(:));
% %     dn = 0.5*dn
% 
% 
% 
% 
% end
% 'counter'
% counter
% 
% 
% 
% if interp
%     'interpolating'
% 
%     if val_c2_prior > val_c1_prior % went with left
% 
%         %ft = fit([-1.5,-0.5,1]',[val_c2, val_c1, val_c1_prior/2]','poly2');
%         ft = polyfit( [-1.5,-0.5,1]',[val_c2, val_c1, val_c1_prior/2]', 2);
%         %mx = -ft.p2/ft.p1/2;
%         mx = -ft(2)/2/ft(1);
%         %if ft.p1<0 && mx > -1.5 && mx < 1
%         if ft(1)<0 && mx > -1.5 && mx < 1
%             p_current = p_c0 + mx*(p_c1_prior - p_c0);
%         end
%         val_c1
%         val_c2
%         val_c1_prior/2
%     else
%         % ft = fit([-1,0.5,1.5]',[val_c2_prior/2, val_c2, val_c1]','poly2');
%         % mx = -ft.p2/ft.p1/2;
%         ft = polyfit([-1,0.5,1.5]',[val_c2_prior/2, val_c2, val_c1]',2);
%         mx = -ft(2)/2/ft(1);
%         %if ft.p1 < 0 && mx > -1 && mx < 1.5
%         if ft(1) < 0 && mx > -1 && mx < 1.5
%             p_current = p_c0 + mx*(p_c1_prior - p_c0);
%         end
% 
% 
%         val_c2_prior/2
%         val_c2
%         val_c1
% 
%     end
% 
% end
% 
% if norm( p_current - [0.053; 0.03] ) < 0.01
%     111
% end


% end