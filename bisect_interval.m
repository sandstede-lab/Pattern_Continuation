
function [p_current,counter, sim_p_history, sim_metric_history, val_c1, val_c2]=bisect_interval(modelpar, featpar, feature_handle, obj_handle, p_c1, p_c2, dn, max_count)
%blablabla

interp = 1;
max_bisect = 2;


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
feature_bp_old = [];
feature_bm_old = [];

p_c0 = p_c1 + p_c2; p_c0 = p_c0/2;

    
    
flag = 0;
counter = 0;
feature_bp = [];feature_bm = [];feature_b0 = [];
val_c2 = 0; val_c1 = 0; 
dists_c2 = []; dists_c1 = []; 
    
  

featpar.N = 1;
modelpar.a = p_c1(1) + dn(1);
modelpar.b = p_c1(2) + dn(2);

counter = 0;
flag = 0;
while (counter < max_count) && (flag == 0)
    
    modelpar.a = p_c1(1) + dn(1);
    modelpar.b = p_c1(2) + dn(2);

    feature_bp{counter+1} = feature_handle(featpar, modelpar);
    if iscell(feature_bp{counter+1})
        feature_bp{counter+1} = feature_bp{counter+1}{1};
    end
    modelpar.a = p_c1(1) - dn(1);%( p_c1(1) + p_c2(1) )/2;
    modelpar.b = p_c1(2) - dn(2);%( p_c1(2) + p_c2(2) )/2;
    
    feature_b0{counter+1} = feature_handle(featpar, modelpar);
    if iscell(feature_b0{counter+1})
        feature_b0{counter+1} = feature_b0{counter+1}{1};
    end
    
    dists_c1(counter+1) = obj_handle(feature_bp{counter+1},feature_b0{counter+1});
    
    modelpar.a = p_c2(1) - dn(1);
    modelpar.b = p_c2(2) - dn(2);
    
    
    feature_bm{counter+1} = feature_evaluation(featpar, modelpar);
    if iscell(feature_bm{counter+1})
        feature_bm{counter+1} = feature_bm{counter+1}{1};
    end
    
    dists_c2(counter+1) = obj_handle(feature_bm{counter+1},feature_b0{counter+1});
    if counter > 1
        flag = ttest(dists_c1(:),dists_c2(:));
    end
    counter = counter + 1;
end
    
    
feature_bm_all = [];
feature_bp_all = [];
feature_b0_all = [];

for ii = 1:length(feature_bm)
        switch featpar.feature
            case 'roundness'
                feature_bm_all{ii} = feature_bm{ii};

                feature_bp_all{ii} = feature_bp{ii};

                feature_b0_all{ii} = feature_b0{ii};
            otherwise
        
                feature_bm_all = [feature_bm_all,feature_bm{ii}];
                
                feature_bp_all = [feature_bp_all,feature_bp{ii}];

                feature_b0_all = [feature_b0_all,feature_b0{ii}];
        end
        
        
        
%                 end
end
        
    
% if max_count > 1
%     flag = ttest(dists_c1(:),dists_c2(:));
% end
        % if counter > 1 || max_count == 1
        
val_c2 = obj_handle(feature_bm_all, feature_b0_all);
val_c1 = obj_handle(feature_bp_all, feature_b0_all);
        
        % end
if featpar.avoid_steady == 1% avoid homogeneous unless we want that
    if feature_bp_all == 0
        val_c1 = -999;
    end
    if feature_bm_all == 0
        val_c2 = -999;
    end
    if feature_b0_all == 0
        val_c1 = -999;
        val_c2 = -999;
    end
end

sim_p_history = [sim_p_history, p_c1'+dn', p_c1'-dn', p_c2'-dn'];

sim_metric_history = [sim_metric_history, mean(cat(2,feature_bp_all)), mean(cat(2,feature_b0_all)), mean(cat(2,feature_bm_all))];
        
        
    
    
    val_c1
    val_c2
    
p_c1_prior = p_c1;
p_c2_prior = p_c2;

dn = dn/2;
            
if val_c2 > val_c1 
%         if steady == 0 
    choice = 'left';
    p_current = p_c2;
    p_prev = p_c1;
    if counter_bisect  < max_bisect
        
        
        
        % if counter_bisect == 1
        %     dn = 0.5*dn;
        % 
        %     p_c2 = p_c2 - dn;
        %     p_c1 = p_c2 + 2*dn;
        %     %p_c1 = p_c2 + dn;
        % end
        
    end

    feature_bp = feature_b0; 
    % p_c2 = p_c2 - dn;
    % p_c1 = p_c2 + 2*dn;
    
    p_c1 = p_c2 + dn;
    p_c2 = p_c2 - dn;
    
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

        % if counter_bisect == 1
        %     dn = 0.5*dn;
        % 
        %     p_c2 = p_c1 - dn;
        %     p_c1 = p_c2 + 2*dn;
        %     %p_c1 = p_c1 + dn;
        % end
    end

    feature_bm = feature_b0; 
    
    % p_c2 = p_c1 - dn;
    % p_c1 = p_c2 + 2*dn;
    p_c2 = p_c1 - dn;
    p_c1 = p_c1 + dn;
    
end
        
val_c1_prior = val_c1;
val_c2_prior = val_c2;


% counter = 0;

flag = 0;
dists_c1 = [];
dists_c2 = [];

%while counter < length(dists_c1)-1 && flag
for counter = 0:length(feature_bm)-1    
    
    modelpar.a = ( p_c1(1) + p_c2(1) )/2;
    modelpar.b = ( p_c1(2) + p_c2(2) )/2;
    
    feature_b0{counter+1} = feature_handle(featpar, modelpar);
    if iscell(feature_b0{counter+1})
        feature_b0{counter+1} = feature_b0{counter+1}{1};
    end
    dists_c1(counter+1) = obj_handle(feature_bp{counter+1},feature_b0{counter+1});
    dists_c2(counter+1) = obj_handle(feature_bm{counter+1},feature_b0{counter+1});
    % if counter > 1
    %     flag = ttest(dists_c1(:),dists_c2(:));
    % end
    %counter = counter + 1;
end


if length(dists_c1) > 1
    flag = ttest(dists_c1(:),dists_c2(:));
end
counter = length(dists_c1) ;

dn = (p_c1 - p_c2)/2;
% counter = 0;
while (counter < max_count) && (flag == 0)
    
    modelpar.a = p_c1(1) + dn(1);
    modelpar.b = p_c1(2) + dn(2);

    feature_bp{counter+1} = feature_handle(featpar, modelpar);
    if iscell(feature_bp{counter+1})
        feature_bp{counter+1} = feature_bp{counter+1}{1};
    end
    modelpar.a = ( p_c1(1) + p_c2(1) )/2;
    modelpar.b = ( p_c1(2) + p_c2(2) )/2;
    
    feature_b0{counter+1} = feature_handle(featpar, modelpar);
    if iscell(feature_b0{counter+1})
        feature_b0{counter+1} = feature_b0{counter+1}{1};
    end
    
    dists_c1(counter+1) = obj_handle(feature_bp{counter+1},feature_b0{counter+1});
    
    modelpar.a = p_c2(1) - dn(1);
    modelpar.b = p_c2(2) - dn(2);
    
    
    feature_bm{counter+1} = feature_evaluation(featpar, modelpar);
    if iscell(feature_bm{counter+1})
        feature_bm{counter+1} = feature_bm{counter+1}{1};
    end
    
    dists_c2(counter+1) = obj_handle(feature_bm{counter+1},feature_b0{counter+1});
    if counter > 1
        flag = ttest(dists_c1(:),dists_c2(:));
    end
    counter = counter + 1;
end




feature_bm_all = [];
feature_bp_all = [];
feature_b0_all = [];

for ii = 1:length(feature_bm)
        switch featpar.feature
            case 'roundness'
                feature_bm_all{ii} = feature_bm{ii};

                feature_bp_all{ii} = feature_bp{ii};

                feature_b0_all{ii} = feature_b0{ii};
            otherwise
        
                feature_bm_all = [feature_bm_all,feature_bm{ii}];
                
                feature_bp_all = [feature_bp_all,feature_bp{ii}];

                feature_b0_all = [feature_b0_all,feature_b0{ii}];
        end
        
        
        
%                 end
end
        
    
% if max_count > 1
%     flag = ttest(dists_c1(:),dists_c2(:));
% end
        % if counter > 1 || max_count == 1
        
val_c2 = obj_handle(feature_bm_all, feature_b0_all);
val_c1 = obj_handle(feature_bp_all, feature_b0_all);
        

sim_p_history = [sim_p_history, p_c1'+dn', p_c1'-dn', p_c2'-dn'];

sim_metric_history = [sim_metric_history, mean(cat(2,feature_bp_all)), mean(cat(2,feature_b0_all)), mean(cat(2,feature_bm_all))];
  
    
    
sim_p_history

'counter'
counter

if featpar.avoid_steady == 1% avoid homogeneous unless we want that
    if feature_bp_all == 0
        val_c1 = -999;
    end
    if feature_bm_all == 0
        val_c2 = -999;
    end
    if feature_b0_all == 0
        val_c1 = -999;
        val_c2 = -999;
    end
end


if interp
    'interpolating'
    
    if val_c2_prior > val_c1_prior % went with left
        
        %ft = fit([-1.5,-0.5,1]',[val_c2, val_c1, val_c1_prior/2]','poly2');
        ft = polyfit( [-1.5,-0.5,1]',[val_c2, val_c1, val_c1_prior/2]', 2);
        %mx = -ft.p2/ft.p1/2;
        mx = -ft(2)/2/ft(1);
        %if ft.p1<0 && mx > -1.5 && mx < 1
        if ft(1)<0 && mx > -1.5 && mx < 1
            p_current = p_c0 + mx*(p_c1_prior - p_c0);
        elseif val_c2 > val_c1
            p_current = p_c2;
        end
        val_c1
        val_c2
        val_c1_prior/2
    else
        % ft = fit([-1,0.5,1.5]',[val_c2_prior/2, val_c2, val_c1]','poly2');
        % mx = -ft.p2/ft.p1/2;
        ft = polyfit([-1,0.5,1.5]',[val_c2_prior/2, val_c2, val_c1]',2);
        mx = -ft(2)/2/ft(1);
        %if ft.p1 < 0 && mx > -1 && mx < 1.5
        if ft(1) < 0 && mx > -1 && mx < 1.5
            p_current = p_c0 + mx*(p_c1_prior - p_c0);
        elseif val_c1 > val_c2
            p_current = p_c1;
        end
        
        
        val_c2_prior/2
        val_c2
        val_c1
        
    end
    
end




