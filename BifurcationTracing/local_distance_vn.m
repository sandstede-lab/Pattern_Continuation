function [out,Ubp,Ubm,feature_bp_all,feature_bm_all] = local_distance_vn(a, b, dn, Rep,alpha,dist,t_max,model,feature_bp_all,feature_bm_all,sets,steady)

if nargin < 6
    dist = 'num';
end
if nargin < 8
    t_max = 100;
end
if nargin < 9
    feature_bp_all = [];
    
end
if nargin < 10
    feature_bm_all = [];
end

if nargin < 11
    sets = 'pos';
end

if nargin < 12
    steady = 0;
end

switch sets
    case 'pos'
        J = 1;
    case 'neg'
        J = 2;

end


feature_bp = cell(1,Rep);
feature_bm = cell(1,Rep);



switch model
    case 'Brusselator'
        Nx = 50;
        [X,Y] = meshgrid( linspace(0,50,Nx), linspace(0,50,Nx));
        D = 50;
        %alpha = 1.2;
    case 'SH'
        [X,Y] = meshgrid( linspace(0,16*pi,128), linspace(0,16*pi,128));
        D = 16*pi;
        %alpha = 0.5;
    case 'GS'
        [X,Y] = meshgrid( linspace(0,2.5,200), linspace(0,2.5,200));
        D = 2.5;
    
    case 'Schnakenberg'
        D = 4;
        [X,Y] = meshgrid( linspace(0,D,80), linspace(0,D,80));
        
end
X = X(:); Y = Y(:);

viz = false;


Ubp = cell(1,Rep);
Ubm = cell(1,Rep);
if strcmp(model,'Brusselator')
    Vbp = cell(1,Rep);
    Vbm = cell(1,Rep);
end


if length(feature_bp_all) == 0
    for j = 1:Rep
        cd ../DataGenerator/Reaction_Diffusion
            switch model
                case 'Brusselator'
                    [Ubp{j},Vbp{j}] = Brusselator(a+dn(1),b+dn(2),viz, Nx,t_max);
                case 'SH'
                    Ubp{j}= SH_2D(a+dn(1),b+dn(2));
                case 'GS'
                    Ubp{j}= Gray_Scott(a+dn(1),b+dn(2));
                case 'Schnakenberg'
                    Ubp{j}= Schnakenberg(a+dn(1),b+dn(2));
            end
        cd ../../BifurcationTracing
        flag =  (max(Ubp{j}(:)) - min(Ubp{j}(:)) > 1E-1);
        if steady == 0 && flag
            idx = find(Ubp{j} > quantile(Ubp{j}(:), 0.7));
            shp_pos = alphaShape(X(idx),Y(idx),alpha);

            idx = find(Ubp{j} < quantile(Ubp{j}(:), 0.3));
            shp_neg = alphaShape(X(idx),Y(idx),alpha);
            figure(28);
            subplot(1,3,1)
            plot(shp_pos)
            xlim([0,D])
            ylim([0,D])
            subplot(1,3,2)
            plot(shp_neg)
            xlim([0,D])
            ylim([0,D])
            subplot(1,3,3)
            imagesc(Ubp{j})
            colorbar()
            title('positive')

            switch dist
                case 'num'
                    feature_bp{j}{1} = numRegions(shp_pos);
                    feature_bp{j}{2} = numRegions(shp_neg);
                case 'area'
                    feature_bp{j}{1} = area(shp_pos, 1:numRegions(shp_pos));
                    p1 = perimeter(shp_pos, 1:numRegions(shp_pos));
                    feature_bp{j}{1} = feature_bp{j}{1}./p1.*min(p1,D);



                    feature_bp{j}{2} = area(shp_neg, 1:numRegions(shp_neg));
                    p2 = perimeter(shp_neg, 1:numRegions(shp_neg));
                    feature_bp{j}{2} = feature_bp{j}{2}./p2.*min(p2,D);
                case 'perimeter'
                    feature_bp{j}{1} = perimeter(shp_pos, 1:numRegions(shp_pos));
                    feature_bp{j}{2} = perimeter(shp_neg, 1:numRegions(shp_neg));

                    feature_bp{j}{1} = min(feature_bp{j}{1} , D);
                    feature_bp{j}{2} = min(feature_bp{j}{2} , D);



                case {'roundness','roundness-ks'}
                    areas = area(shp_pos, 1:numRegions(shp_pos));
                    perimeters = perimeter(shp_pos, 1:numRegions(shp_pos));   
                    feature_bp{j}{1} = 4*pi*areas./(perimeters.^2);


                    areas = area(shp_neg, 1:numRegions(shp_neg));
                    perimeters = perimeter(shp_neg, 1:numRegions(shp_neg));   
                    feature_bp{j}{2} = 4*pi*areas./(perimeters.^2);

                case 'PCA-roundness'
                    feature_bp{j}{1} = zeros( 1,numRegions(shp_pos) );
                    P = shp_pos.Points;
                    for i = 1:numRegions(shp_pos)
                        out = boundaryFacets(shp_pos,i); out = unique(out(:));
                        samps = P(out,:);

                        out = eigs( (samps - mean(samps))*(samps - mean(samps))');
                        feature_bp{j}{1}(i) = out(2)/out(1);
                    end

                    feature_bp{j}{2} = zeros( 1,numRegions(shp_neg) );
                    P = shp_neg.Points;
                    for i = 1:numRegions(shp_neg)
                        out = boundaryFacets(shp_neg,i); out = unique(out(:));
                        samps = P(out,:);

                        out = eigs( (samps - mean(samps))*(samps - mean(samps))');
                        feature_bp{j}{2}(i) = out(2)/out(1);
                    end



            end
        else
            feature_bp{j}{1} = flag;feature_bp{j}{2} = flag;
            
        end
    end
end
% fprintf('Ubp done\n')
if length(feature_bm_all) == 0
    for j = 1:Rep
        cd ../DataGenerator/Reaction_Diffusion
            switch model
                case 'Brusselator'
                    [Ubm{j},Vbm{j}] = Brusselator(a-dn(1),b-dn(2),viz, Nx,t_max);%, Ubm{j},Vbm{j});
                case 'SH'
                    Ubm{j}= SH_2D(a-dn(1),b-dn(2));
                case 'GS'
                    Ubm{j}= Gray_Scott(a-dn(1),b-dn(2));
                case 'Schnakenberg'
                    Ubm{j}= Schnakenberg(a-dn(1),b-dn(2));
            end
        cd ../../BifurcationTracing
        flag =  (max(Ubm{j}(:)) - min(Ubm{j}(:)) > 1E-1);
        if steady==0 && flag
    
            idx = find(Ubm{j} > quantile(Ubm{j}(:), 0.7));
            shp_pos = alphaShape(X(idx),Y(idx),alpha);

            idx = find(Ubm{j} < quantile(Ubm{j}(:), 0.3));
            shp_neg = alphaShape(X(idx),Y(idx),alpha);
            
            figure(28);
            subplot(1,3,1)
            plot(shp_pos)
            xlim([0,D])
            ylim([0,D])
            subplot(1,3,2)
            plot(shp_neg)
            xlim([0,D])
            ylim([0,D])
            subplot(1,3,3)
            imagesc(Ubm{j})
            colorbar()
            title('negative')
            
            
            switch dist
                case 'num'
                    feature_bm{j}{1} = numRegions(shp_pos);
                    feature_bm{j}{2} = numRegions(shp_neg);
                case 'area'
                    feature_bm{j}{1} = area(shp_pos, 1:numRegions(shp_pos));
                    p1 = perimeter(shp_pos, 1:numRegions(shp_pos));
                    feature_bm{j}{1} = feature_bm{j}{1}./p1.*min(p1,D);



                    feature_bm{j}{2} = area(shp_neg, 1:numRegions(shp_neg));
                    p2 = perimeter(shp_neg, 1:numRegions(shp_neg));
                    feature_bm{j}{2} = feature_bm{j}{2}./p2.*min(p2,D);

                case 'perimeter'
                    feature_bm{j}{1} = perimeter(shp_pos, 1:numRegions(shp_pos));
                    feature_bm{j}{2} = perimeter(shp_neg, 1:numRegions(shp_neg));

                    feature_bm{j}{1} = min(feature_bm{j}{1} , D);
                    feature_bm{j}{2} = min(feature_bm{j}{2} , D);

                case {'roundness','roundness-ks'}
                    areas = area(shp_pos, 1:numRegions(shp_pos));
                    perimeters = perimeter(shp_pos, 1:numRegions(shp_pos));   
                    feature_bm{j}{1} = 4*pi*areas./(perimeters.^2);


                    areas = area(shp_neg, 1:numRegions(shp_neg));
                    perimeters = perimeter(shp_neg, 1:numRegions(shp_neg));   
                    feature_bm{j}{2} = 4*pi*areas./(perimeters.^2);

                case 'PCA-roundness'
                    feature_bm{j}{1} = zeros( 1,numRegions(shp_pos) );
                    P = shp_pos.Points;
                    for i = 1:numRegions(shp_pos)
                        out = boundaryFacets(shp_pos,i); out = unique(out(:));
                        samps = P(out,:);

                        out = eigs( (samps - mean(samps))*(samps - mean(samps))');
                        feature_bm{j}{1}(i) = out(2)/out(1);
                    end

                    feature_bm{j}{2} = zeros( 1,numRegions(shp_neg) );
                    P = shp_neg.Points;
                    for i = 1:numRegions(shp_neg)
                        out = boundaryFacets(shp_neg,i); out = unique(out(:));
                        samps = P(out,:);

                        out = eigs( (samps - mean(samps))*(samps - mean(samps))');
                        feature_bm{j}{2}(i) = out(2)/out(1);
                    end

            end
        else
             feature_bm{j}{1} = flag;feature_bm{j}{2} = flag;
            
        end
    end
end
% fprintf('Ubm done\n')
% feature_bp{1}
% feature_bm{1}

out = 0;
switch dist
    case 'num'

        if isempty(feature_bp_all)
            feature_bp_all = cell(1,length(J));
            for i = 1:Rep
                for j = 1:length(J)
                    feature_bp_all{j} = [feature_bp_all{j},feature_bp{i}{J(j)}];
                end
               
               
            end
        end
        if isempty(feature_bm_all)
            feature_bm_all = cell(1,length(J));
            for i = 1:Rep
                for j = 1:length(J)
                    feature_bm_all{j} = [feature_bm_all{j},feature_bm{i}{J(j)}];
                end
            end
        end

        
        for j = 1:length(J)
            out = out + abs(sum(feature_bp_all{j}) - sum(feature_bm_all{j}))/Rep;
        end

       
                
    otherwise
%         feature_all = [];
        if isempty(feature_bp_all)
            feature_bp_all = cell(1,length(J));
            for i = 1:Rep
    
                for j = 1:length(J)
                    feature_bp_all{j} = [feature_bp_all{j},feature_bp{i}{J(j)}];
                end
               
               
            end
        end
        if isempty(feature_bm_all)
            feature_bm_all = cell(1,length(J));
            for i = 1:Rep
                for j = 1:length(J)
                    feature_bm_all{j} = [feature_bm_all{j},feature_bm{i}{J(j)}];
                end
            end
        end
           
        
        for j = 1:length(J)
            switch dist
                case 'roundness'
                    out = out + ws_distance( feature_bm_all{j},feature_bp_all{j}, 2 );
                case 'roundness-ks'
                    [~,~,ks] = kstest2( feature_bm_all{j},feature_bp_all{j} );
                    out = out + ks;
                case 'PCA-roundness'
                    out = out + ws_distance( feature_bm_all{j},feature_bp_all{j}, 2 );
            end
        end

        
end


