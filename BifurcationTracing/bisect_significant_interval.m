function [p_current,counter, sim_p_history, sim_metric_history] = bisect_significant_interval(model,p_c1, p_c2, dn, max_bisect, alpha,dist,t_max,sets,steady, interp,mode)

switch model
    case 'Brusselator'
        max_count = 10;   
        alpha = 1;
    case 'SH'
        alpha = 0.3+0.05;
        max_count = 5;
    case 'GS'
        alpha = 0.01;
        max_count = 5;
    case 'Oregonator'
        alpha = 1.5;
        max_count = 5;
    case 'Schnakenberg'
        alpha = 0.1;
%         alpha = 0.08;
        max_count = 10;
end

if steady
    max_count = 1;
end

switch sets
    case 'pos'
        J = 1;
    case 'neg'
        J = 2;
    case 'both'
        J = [1,2];
end
if nargin < 11
    interp = 0;
end


counter_bisect = 0;
flag = 0;
counter = 0;
feature_bp_all = [];
feature_bm_all = []; 
feature_b0_all = [];
sim_p_history = [];
sim_metric_history = [];

counter = 0;

feature_bp = [];feature_bm = [];feature_b0 = [];
val_c2 = 0; val_c1 = 0; 
dists_c2 = []; dists_c1 = []; 
feature_bp_old = []
feature_bm_old = [];

p_c0 = p_c1 + p_c2; p_c0 = p_c0/2;
while counter_bisect < max_bisect
    
    
    flag = 0;
    counter = 0;
    feature_bp = [];feature_bm = [];feature_b0 = [];
    val_c2 = 0; val_c1 = 0; 
    dists_c2 = []; dists_c1 = []; 
    
  

    while flag == 0
       
        
        if counter < max_count
            if counter+1 > length(feature_bp_old) || counter_bisect == 0
                [dists_c1(counter+1),~,~,feature_bp{counter+1},feature_b0{counter+1}] = local_distance_vn(p_c1(1), p_c1(2), dn, 1,alpha,dist,t_max,model,[],[],sets,steady);
            
                [dists_c2(counter+1),~,~,~,feature_bm{counter+1}] = local_distance_vn(p_c2(1), p_c2(2), dn, 1,alpha,dist,t_max,model,feature_b0{counter+1},[],sets,steady);
            else
                feature_bp{counter+1} = feature_bp_old{counter+1};
                [dists_c1(counter+1),~,~,~,feature_b0{counter+1}] = local_distance_vn(p_c1(1), p_c1(2), dn, 1,alpha,dist,t_max,model,feature_bp{counter+1},[],sets,steady);
                feature_bm{counter+1} = feature_bm_old{counter+1};
                [dists_c2(counter+1),~,~,~,~] = local_distance_vn(p_c2(1), p_c2(2), dn, 1,alpha,dist,t_max,model,feature_b0{counter+1},feature_bm{counter+1},sets,steady);
            end
%         end
            counter = counter + 1;
%         else
%             counter = max_count;
        
        end
        
        val_c1 = 0; val_c2=0; 
        
        
        feature_bm_all = cell(1,length(J)); 
        feature_bp_all = cell(1,length(J)); 
        feature_b0_all = cell(1,length(J)); 
        for j = 1:length(J)
            for ii = 1:counter
                if strcmp(mode, 'dist')
                
                    feature_bm_all{j} = [feature_bm_all{j},feature_bm{ii}{j}];

                    feature_bp_all{j} = [feature_bp_all{j},feature_bp{ii}{j}];

                    feature_b0_all{j} = [feature_b0_all{j},feature_b0{ii}{j}];
                else
                    
                    feature_bm_all{j} = [feature_bm_all{j},feature_bm{ii}{j}];

                    feature_bp_all{j} = [feature_bp_all{j},feature_bp{ii}{j}];

                    feature_b0_all{j} = [feature_b0_all{j},feature_b0{ii}{j}];
                    
                    
                    feature_bm_all_cell{j}{ii} = feature_bm{ii}{j};

                    feature_bp_all_cell{j}{ii} = feature_bp{ii}{j};

                    feature_b0_all_cell{j}{ii} = feature_b0{ii}{j};
                end
            end
        end
%         feature_bm_all{1}
%         feature_bp_all{1}
%         feature_b0_all{1}
        
%         [~,idx] = sort([mean(dists_c1(:)),mean(dists_c2(:))]);
        
        if counter > 1
            flag = ttest(dists_c1(:),dists_c2(:));
        end
        if counter > 1 || max_count == 1
            for j = 1:length(J)
                if strcmp( mode, 'dist')
                    val_c2 = val_c2 + feature_distance(cat(2,feature_bm_all{j}), cat(2,feature_b0_all{j}), dist, 'dist'); 
                    val_c1 = val_c1 + feature_distance(cat(2,feature_b0_all{j}), cat(2,feature_bp_all{j}), dist, 'dist'); 
                else
                    if strcmp(dist,'roundness')
                        val_c2 = val_c2 + feature_distance(feature_bm_all_cell{j}, feature_b0_all_cell{j}, dist, mode);
                        val_c1 = val_c1 + feature_distance(feature_b0_all_cell{j}, feature_bp_all_cell{j}, dist, mode);
                    else
                        val_c2 = val_c2 + feature_distance(feature_bm_all{j}, feature_b0_all{j}, dist, mode);
                        val_c1 = val_c1 + feature_distance(feature_b0_all{j}, feature_bp_all{j}, dist, mode);
                    end
                end
                
            end
        end
        if steady == 0 % avoid homogeneous unless we want that
            if feature_bp{1}{1}(1) == 0
                val_c1 = -999;
            end
            if feature_bm{1}{1}(1) == 0
                val_c2 = -999;
            end
        end
        
        sim_p_history = [sim_p_history, p_c1'+dn', p_c1'-dn', p_c2'-dn'];
        
        sim_metric_history = [sim_metric_history, mean(cat(2,feature_bp_all{1})), mean(cat(2,feature_b0_all{1})), mean(cat(2,feature_bm_all{1}))];
        
        if counter >= max_count
            flag = 1;
        end
    end
    
    
    
    
%     if counter_bisect + 1 < max_bisect
%         p_c2_prior = p_c2; 
%         p_c1_prior = p_c1;
        counter_bisect = counter_bisect + 1
        if counter_bisect == 1
            val_c2_prior = val_c2;
            val_c1_prior = val_c1;
            p_c2_prior = p_c2; 
            p_c1_prior = p_c1;
            
        end
            
        if val_c2 > val_c1 
    %         if steady == 0 
                choice = 'left';
                p_current = p_c2;
                p_prev = p_c1;
                if counter_bisect  < max_bisect
                    
                    
                    
                    if counter_bisect == 1
                        dn = 0.5*dn;
                        
                        p_c2 = p_c2 - dn;
                        p_c1 = p_c2 + 2*dn;
                    end
                    
                end

                feature_bp_old = feature_b0; 
                feature_bm_old = feature_bm;

    %         else
    %             
    %         end


        else
            choice = 'right';
            p_current = p_c1;
            p_prev = p_c2;
            if counter_bisect  < max_bisect
                
                val_c2_prior = val_c2;
                val_c1_prior = val_c1;
%                 p_c2_prior = p_c2; 
%                 p_c1_prior = p_c1;
                if counter_bisect == 1
                    dn = 0.5*dn;
%                     p_c2_prior = p_c2; 
%                     p_c1_prior = p_c1;
                    p_c2 = p_c1 - dn;
                    p_c1 = p_c2 + 2*dn;
                end
            end

            feature_bm_old = feature_b0; 
            feature_bp_old = feature_bp;
        end
        
%     end
    
%     flag = ttest(dists_c1(:),dists_c2(:));
%     dn = 0.5*dn
    
    
    
    
end
'counter'
counter
if interp
    'interpolating'
    
    if val_c2_prior > val_c1_prior % went with left
        
        ft = fit([-1.5,-0.5,1]',[val_c2, val_c1, val_c1_prior/2]','poly2');
        mx = -ft.p2/ft.p1/2;
        if ft.p1<0 && mx > -1.5 && mx < 1
            p_current = p_c0 + mx*(p_c1_prior - p_c0);
        end
        
    else
        ft = fit([-1,0.5,1.5]',[val_c2_prior/2, val_c2, val_c1]','poly2');
        mx = -ft.p2/ft.p1/2;
        if ft.p1 < 0 && mx > -1 && mx < 1.5
            p_current = p_c0 + mx*(p_c1_prior - p_c0);
        end
        
    end
    
end


    
    
    
    
    
    
    
    
    
    
