function [out,Ubp,Ubm,feature_bp_all,feature_bm_all] = local_distance_vn(a, b, dn, Rep,alpha,dist,t_max,model,feature_bp_all,feature_bm_all,sets,steady)

if nargin < 5
    dist = 'num';
end
if nargin < 6
    t_max = 100;
end
if nargin < 7
    model = 'Brusselator';
end
if nargin < 8
    feature_bp_all = [];
    
end
if nargin < 9
    feature_bm_all = [];
end

if nargin < 10
    sets = 'pos';
end

if nargin < 11
    steady = 0;
end




feature_bp = cell(1,Rep);
feature_bm = cell(1,Rep);



switch model
    case 'Brusselator'
        Nx = 50;
        [X,Y] = meshgrid( linspace(0,50,Nx), linspace(0,50,Nx));
        D = 50;
        %alpha = 1.2;
        X = X(:); Y = Y(:);
    case 'SH'
        [X,Y] = meshgrid( linspace(0,16*pi,128), linspace(0,16*pi,128));
        D = 16*pi;
        X = X(:); Y = Y(:);
%     case 'SH_1D'
        
    case 'GS'
        [X,Y] = meshgrid( linspace(0,2.5,200), linspace(0,2.5,200));
        D = 2.5;
        X = X(:); Y = Y(:);
    case 'Schnakenberg'
        D = 4;
        [X,Y] = meshgrid( linspace(0,D,80), linspace(0,D,80));
        X = X(:); Y = Y(:);
end
% X = X(:); Y = Y(:);

viz = false;


Ubp = cell(1,Rep);
Ubm = cell(1,Rep);
if strcmp(model,'Brusselator')
    Vbp = cell(1,Rep);
    Vbm = cell(1,Rep);
end


if isempty(feature_bp_all)
    for j = 1:Rep
%         cd ../DataGenerator/Reaction_Diffusion
            switch model
                case 'Brusselator'
                    [Ubp{j},Vbp{j}] = Brusselator(a+dn(1),b+dn(2),viz, Nx,t_max);
                case 'SH'
                    Ubp{j}= SH_2D(a+dn(1),b+dn(2));
                case 'SH_1D'
                    Ubp{j}= SH_1D(a+dn(1),b+dn(2));
                case 'GS'
                    Ubp{j}= Gray_Scott(a+dn(1),b+dn(2));
                case 'Schnakenberg'
                    Ubp{j}= Schnakenberg(a+dn(1),b+dn(2));
            end
%         cd ../../BifurcationTracing
        flag =  (max(Ubp{j}(:)) - min(Ubp{j}(:)) > 1E-1);
        if steady == 0 && flag 
            if min( size(Ubp{j}) ) > 1
                switch sets
                    case 'pos'
                        idx = find(Ubp{j} > quantile(Ubp{j}(:), 0.7));
                    case 'neg'
                        idx = find(Ubp{j} < quantile(Ubp{j}(:), 0.3));
                end
                shp = alphaShape(X(idx),Y(idx),alpha);

                
                figure(28);
                subplot(1,2,1)
                plot(shp)
                xlim([0,D])
                ylim([0,D])
                
                subplot(1,2,2)
                imagesc(Ubp{j})
                colorbar()
                title('positive')
            
                
                
            end

            switch dist
                case 'num'
                    if min( size(Ubp{j}) ) > 1
                        feature_bp{j} = numRegions(shp);
                        
                    else
                        
                        idx = find(Ubp{j}(:) > 0.5*(max(Ubp{j}(:)) ));
                        feature_bp{j} = sum( diff(idx)>1 );
                         
                    end
               
                case {'roundness','roundness-bag'}
                    areas = area(shp, 1:numRegions(shp));
                    perimeters = perimeter(shp, 1:numRegions(shp));   
                    feature_bp{j} = 4*pi*areas./(perimeters.^2);


                


            end
        
            
        else
            feature_bp{j} = flag;
            
        end
    end
end
% fprintf('Ubp done\n')
if isempty(feature_bm_all)
    for j = 1:Rep
%         cd ../DataGenerator/Reaction_Diffusion
            switch model
                case 'Brusselator'
                    [Ubm{j},Vbm{j}] = Brusselator(a-dn(1),b-dn(2),viz, Nx,t_max);%, Ubm{j},Vbm{j});
                case 'SH'
                    Ubm{j}= SH_2D(a-dn(1),b-dn(2));
                case 'SH_1D'
                    Ubm{j}= SH_1D(a-dn(1),b-dn(2));
                case 'GS'
                    Ubm{j}= Gray_Scott(a-dn(1),b-dn(2));
                case 'Schnakenberg'
                    Ubm{j}= Schnakenberg(a-dn(1),b-dn(2));
            end
%         cd ../../BifurcationTracing
        flag =  (max(Ubm{j}(:)) - min(Ubm{j}(:)) > 1E-1);
        if steady==0 && flag
            if min( size(Ubm{j}) ) > 1
                switch sets
                    case 'pos'
                        idx = find(Ubm{j} > quantile(Ubm{j}(:), 0.7));
                    case 'neg'
                        idx = find(Ubm{j} < quantile(Ubm{j}(:), 0.3));
                end
                shp = alphaShape(X(idx),Y(idx),alpha);

                
                figure(28);
                subplot(1,2,1)
                plot(shp)
                xlim([0,D])
                ylim([0,D])
                
                subplot(1,2,2)
                imagesc(Ubm{j})
                colorbar()
                title('negative')
            end
            
            
            switch dist
                case 'num'
                    if min( size(Ubm{j}) ) > 1
                        feature_bm{j} = numRegions(shp);
                        
                    else
                        idx = find(Ubm{j}(:) > 0.5*(max(Ubm{j}(:)) ));
                        feature_bm{j} = sum( diff(idx)>1 );
                        
                    end
               
                case {'roundness','roundness-bag'}
                    areas = area(shp, 1:numRegions(shp));
                    perimeters = perimeter(shp, 1:numRegions(shp));   
                    feature_bm{j} = 4*pi*areas./(perimeters.^2);



                
            end
        else
             feature_bm{j} = flag;
            
        end
    end
end
% fprintf('Ubm done\n')
% feature_bp{1}
% feature_bm{1}


switch dist
    case 'num'

        if isempty(feature_bp_all)
            feature_bp_all = [];
            for i = 1:Rep
                
                feature_bp_all = [feature_bp_all,feature_bp{i}];
                
               
               
            end
        end
        if isempty(feature_bm_all)
            feature_bm_all = [];
            for i = 1:Rep
                
                feature_bm_all = [feature_bm_all,feature_bm{i}];
                
               
               
            end
        end

        
        out = abs(sum(feature_bp_all) - sum(feature_bm_all))/Rep;
        

       
                
    otherwise
%         feature_all = [];
        if isempty(feature_bp_all)
            feature_bp_all = [];
            for i = 1:Rep
    
                
                    feature_bp_all = [feature_bp_all,feature_bp{i}];
               
               
               
            end
        end
        if isempty(feature_bm_all)
            feature_bm_all = [];
            for i = 1:Rep
                
                    feature_bm_all = [feature_bm_all,feature_bm{i}];
                
            end
        end
           
        
        out = ws_distance( feature_bm_all,feature_bp_all, 2 );
            
        

        
end


