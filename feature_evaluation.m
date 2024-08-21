function features = feature_evaluation(featpar, modelpar)

% simulate the patterns to be used for statistics
pattern = cell(1,featpar.N);
parfor i = 1:featpar.N
    pattern{i} = model(modelpar);
end

switch featpar.feature
    case {'num','roundness-bag','steady','retract','meander','drift','turbulence','area'}
        features = [];
    otherwise
        features = cell(1,featpar.N);
end

% compute statistics based on patterns

for i = 1:featpar.N

    if strcmp( modelpar.sets, 'tip_points' )
        if ~isempty(pattern{i})
            X = pattern{i}(:,2);
            Y = pattern{i}(:,3);
        
        end

    else

        if ~isempty(pattern{i})
            X = pattern{i}(:,1);
        end
        
        if modelpar.xdim == 2 && ~isempty(pattern{i})
            Y = pattern{i}(:,2);
        
            shp = alphaShape(X,Y,featpar.alpha);
        
        end
    end
    
    switch featpar.feature
        case {'steady','retract'}
            if ~isempty(pattern{i})
                features = [features, 1];
            else
                features = [features, 0];
            end
        case 'num'
            if ~isempty(pattern{i})
                switch modelpar.xdim
                    case 2
                        features = [features, numRegions(shp)];
                    
                    case 1
                    
                        features = [features, sum( diff(X)>min(diff(X)) )];
                     
                end
            else
                features = [features,0];
            end
       
        case 'roundness-bag'
            if ~isempty(pattern{i})
    
                areas = area(shp, 1:numRegions(shp));
                perimeters = perimeter(shp, 1:numRegions(shp));   
                features = [features, 4*pi*areas./(perimeters.^2)];
            end
        case 'roundness'
            if ~isempty(pattern{i})
                areas = area(shp, 1:numRegions(shp));
                perimeters = perimeter(shp, 1:numRegions(shp)); 
                features{i} = 4*pi*areas./(perimeters.^2);
            end

        

        case 'meander'
            if size(pattern{i},1) < 100
                features = 0;
            else

                alpha0 = max( Y )-min(Y);
                alpha0 = max( max( X)-min(X),alpha0);
                alpha = alpha0/2.5;
           
                shp_pos = alphaShape(X,Y,alpha);%,alpha,'RegionThreshold',1);
                
                shp_area = sum(area(shp_pos,1:numRegions(shp_pos)));
                shp_rdns = shp_area/(max( Y )-min(Y));
                shp_rdns = shp_rdns/(max( X)-min(X));
%                         
                features = [features, tanh(4*shp_rdns)];
                
                figure(11);
                plot(shp_pos)
                title( features )
            end

        case 'drift'
            % if size(pattern{i},1) < 100
            %       features = 0;
            % end
            % pattern{i}
            tip_pts = pattern{i}(:,2:3);
            out = eigs( ( tip_pts- mean(tip_pts))'*(tip_pts - mean(tip_pts)));
            l = out(2)/out(1);
            
            if size(pattern{i},1) < 1% || l > 0.99
                features = 0;
                pattern{i} = [0,0,0];
                tip_pts = [0,0];
            else
                while l > 0.6
                    tip_pts = tip_pts(1:round(size(tip_pts,1)*0.75),:);
                    out = eigs( ( tip_pts- mean(tip_pts))'*(tip_pts - mean(tip_pts)));
                    l = out(2)/out(1);
                end
            
                alpha = max( tip_pts(:,2))-min(tip_pts(:,2));
                alpha = max( max( tip_pts(:,1))-min(tip_pts(:,1)),alpha);
                
    %                         while abs(max( tip_pts(:,2))-min(tip_pts(:,2)) - max( tip_pts(:,1))+min(tip_pts(:,1))) < 0.01*alpha
    %                             tip_pts = tip_pts(1:round(size(tip_pts,1)*0.75),:);
    %                             alpha = max( tip_pts(:,2))-min(tip_pts(:,2));
    %                             alpha = max( max( tip_pts(:,1))-min(tip_pts(:,1)),alpha);
    %                         
    %                         end
                
                shp_pos = alphaShape(tip_pts(:,1),tip_pts(:,2),10*alpha,'RegionThreshold',1);
                if numRegions(shp_pos) < 1 && size(pattern{i},1)>1
                    
                    shp_pos = alphaShape(tip_pts(:,1),tip_pts(:,2));%,alpha,'RegionThreshold',1);
                    pc = criticalAlpha(shp_pos,'one-region');
                    shp_pos.Alpha = pc;
                end
                [~,X] = boundaryFacets(shp_pos);
                dists = sqrt(sum( (X-X([2:end,1],:)).^2,2) );
                idx = find(dists == max(dists));
                if idx < size(X,1)
                    p0 = X(idx,:); p1=X(idx+1,:);
                else
                    p0 = X(idx,:); p1=X(1,:);
                end
                if norm(p0-tip_pts(1,:)) > norm(p1-tip_pts(1,:))
                    p_temp = p1; p1 = p0; p0 = p_temp;
                end
                p_dir = p1 - p0;p_dir = [p_dir(2),-p_dir(1)];
            

            
                features = mean(  sign((tip_pts - p0)*p_dir'));
            end
            tip_pts


        case 'turbulence'
            tip_times = pattern{i}(:,1);
            % if length(tip_times)>1000
            %     tip_times = tip_times(end-1000:end);
            %     pattern{i} = pattern{i}(end-1000:end,:);
            % end
            if length(unique(tip_times)) == length(tip_times)
                features = [features,0];
            else
                time = mode(tip_times);
                idx = find( tip_times == time);
                features = 1*( max( reshape(pdist( pattern{i}(idx,2:3)' ),1,[]) )>10) ;
            end
            
        case 'area'
            % this is only for Rossler
            

            if ~isempty(pattern{i})
    
                areas = area(shp, 1:numRegions(shp));
                
                features = [features, sum(areas)];
            end
    
    
    end
end
'features'
features
