function [out,Ubp,Ubm,feature_bp,feature_bm] = local_distance_vn_spiral_c_alpha(a, b, dn, Rep,alpha,dist,t_max,model,feature_bp,feature_bm,sets,steady)

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
    sets = 'both';
end

if nargin < 12
    steady = 0;
end







Ubp = cell(1,Rep);
Ubm = cell(1,Rep);



if length(feature_bp) == 0
    for j = 1:Rep
        
        switch model
            case 'Barkley'
                anew = num2str(a+dn(1),'%1.6f');
                bnew = num2str(b+dn(2),'%1.6f');
                
                
                length(anew)
                length(bnew)
                
                fid = fopen('../DataGenerator/Spiral_Wave/Barkley/task.txt');
                C=textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                C{1}{1} = anew; C{1}{2} = bnew;
                
                writecell(C{1},'../DataGenerator/Spiral_Wave/Barkley/task.txt','QuoteStrings',0)
%                 if b > 0.05
%                     C{1}{11}(1:5) = num2str(100000);
%                 end
                cd ../DataGenerator/Spiral_Wave/Barkley
                system('rm tip.txt')
                system('touch tip.txt')
                system('./ezspiral')
                
                cd ../../../BifurcationTracing/
                tipping_pts{j} = readtable('../DataGenerator/Spiral_Wave/Barkley/tip.txt');
                u = readtable('../DataGenerator/Spiral_Wave/Barkley/fc.txt');
                u = u.Var2;
                u = reshape(u(5:end),501,501);
            case 'Bar-Eiswirth'
                if a+dn(1) < 0
                    dn(1) = 1E-14 - a;
                end
                anew = num2str(1/(a+dn(1)),'%1.6f');
                bnew = num2str(b+dn(2),'%1.6f');
                
                length(anew)
                length(bnew)
                
                fid = fopen('../DataGenerator/Spiral_Wave/Bar-Eiswirth/task.txt');
                C=textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                C{1}{2} = bnew; C{1}{3} = anew;
                writecell(C{1},'../DataGenerator/Spiral_Wave/Bar-Eiswirth/task.txt','QuoteStrings',0)
                
                cd ../DataGenerator/Spiral_Wave/Bar-Eiswirth
                system('rm tip.txt')
                system('touch tip.txt')
                system('./ezspiral')
                
                cd ../../../BifurcationTracing/
                tipping_pts{j} = readtable('../DataGenerator/Spiral_Wave/Bar-Eiswirth/tip.txt');
                
                u = readtable('../DataGenerator/Spiral_Wave/Bar-Eiswirth/fc.txt');
                u = u.Var2;
                u = reshape(u(5:end),501,501);
            case 'Rossler'
                
                
                cnew = num2str(a+dn(1),'%1.6f');
                anew = num2str(b+dn(2),'%1.6f');
                
                fid = fopen('../DataGenerator/Spiral_Wave/Rossler/task.txt');
                C=textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                C{1}{3} = cnew; C{1}{1} = anew;
                writecell(C{1},'../DataGenerator/Spiral_Wave/Rossler/task.txt','QuoteStrings',0)

                cd ../DataGenerator/Spiral_Wave/Rossler

                system('./ezspiral')
                cd ../../../BifurcationTracing/

                u = readtable('../DataGenerator/Spiral_Wave/Rossler/fc.txt');

                
                u = u.Var2; 
                u = reshape(u(5:end),526,526);
%                 cd ../
                tipping_pts{j} = [];
                
        end
                if size(tipping_pts{j},1) == 0
                    tipping_pts{j} = [0,0];
                else
                	if strcmp(dist,'turbulence')
                        tipping_times = [tipping_pts{j}.Var1];
                    
                        
                    end
                    tipping_pts{j} = [tipping_pts{j}.Var2,tipping_pts{j}.Var3];
                end
                switch model
                    case 'Barkley'
                        if size(tipping_pts{j},1) > 1000 || size(tipping_pts{j},1)==1
                            tipping_pts{j} = tipping_pts{j}(1000:end,:);

                        end
                        idx = find( (tipping_pts{j}(:,1)>98)+(tipping_pts{j}(:,1)<2) +...
                            + (tipping_pts{j}(:,2)>98) + (tipping_pts{j}(:,2)<2)  );
                        
                    case 'Bar-Eiswirth'
                        if strcmp(dist,'turbulence') == 0
                            if size(tipping_pts{j},1) > 2000 || size(tipping_pts{j},1)==1
                                tipping_pts{j} = tipping_pts{j}(2000:end,:);

                            end
                            idx = find( (tipping_pts{j}(:,1)>48)+(tipping_pts{j}(:,1)<2) +...
                                + (tipping_pts{j}(:,2)>48) + (tipping_pts{j}(:,2)<2)  );
                        else
                            idx = [];
                        end
                    case 'Rossler'
                        tipping_pts{j} = [];
                        idx = [];
                        
                end
                if length(idx)
                    tipping_pts{j} = tipping_pts{j}(1:idx(1),:);
                end
                if size(tipping_pts{j},1) && size(tipping_pts{j},2)>1
%                     idx = find( pdist2(tipping_pts{j} , tipping_pts{j}(1,:)) < 0.5 );
%                     idx = idx( find(idx>200));
%                     if length(idx)
%                         tipping_pts{j} = tipping_pts{j}(1:idx(end),:);
%                     end
                    figure(10);
                    scatter(tipping_pts{j}(:,1),tipping_pts{j}(:,2),[],1:length(tipping_pts{j}(:,1)))
                    hold on,
                    scatter(tipping_pts{j}([1,end],1),tipping_pts{j}([1,end],2),'k*')
                    hold off
                end
       
        
            switch dist
                case 'steady'
                   feature_bp{j} = min(1,size(tipping_pts{j},1)/1000);%2*steady*(1 - floor(size(tipping_pts{j},1)/1000));
                   feature_bp{j} = size(tipping_pts{j},1) < 1000;
                   feature_bp{j} = max(u(:))-min(u(:)) < 1E-2;
                     
                case 'ring_fit'
                   
                     if size(tipping_pts{j},1) < 100
                          feature_bp{j} = 0;
                     else
%                         shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2));%,alpha,'RegionThreshold',1);
%                         pc = criticalAlpha(shp_pos,'one-region');
%                         shp_pos.Alpha = pc;
                        alpha0 = max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2));
                        alpha0 = max( max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)),alpha0);
                        alpha = alpha0/2.5;
%                         shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha,'RegionThreshold',1,'HoleThreshold',1);
%                         while numRegions(shp_pos) < 1
%                             alpha = alpha*0.9;
%                             shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha);%,alpha,'RegionThreshold',1);
%                             
%                             
%                         end
%                         shp_area = area(shp_pos,1);
%                         shp_peri = perimeter(shp_pos,1)/2/pi;
% %                         shp_rdns = shp_area/(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)));
% %                         shp_rdns = shp_rdns/(max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)));
%                         shp_rdns = shp_area/pi/shp_peri/shp_peri;
                        
                        shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha);%,alpha,'RegionThreshold',1);
                        
                        shp_area = sum(area(shp_pos,1:numRegions(shp_pos)));
                        shp_rdns = shp_area/(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)));
                        shp_rdns = shp_rdns/(max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)));
%                         

                        feature_bp{j} = tanh(4*shp_rdns);
                        figure(11);
                        plot(shp_pos)
                        title( feature_bp{j} )
                     end
                    
                 case 'line_fit'
                        if size(tipping_pts{j},1) < 100
                              feature_bp{j} = 0;
                        else
                        out = eigs( ( tipping_pts{j}- mean(tipping_pts{j}))'*(tipping_pts{j} - mean(tipping_pts{j})));
                        l = out(2)/out(1);
                        l
                        while l > 0.6
                            tipping_pts{j} = tipping_pts{j}(1:round(size(tipping_pts{j},1)*0.75),:);
                            out = eigs( ( tipping_pts{j}- mean(tipping_pts{j}))'*(tipping_pts{j} - mean(tipping_pts{j})));
                            l = out(2)/out(1);
                        end
                        
                        alpha = max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2));
                        alpha = max( max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)),alpha);
                        
%                         while abs(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)) - max( tipping_pts{j}(:,1))+min(tipping_pts{j}(:,1))) < 0.01*alpha
%                             tipping_pts{j} = tipping_pts{j}(1:round(size(tipping_pts{j},1)*0.75),:);
%                             alpha = max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2));
%                             alpha = max( max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)),alpha);
%                         
%                         end
                        
                        shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),10*alpha,'RegionThreshold',1);
                        if numRegions(shp_pos) < 1
                            
                            shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2));%,alpha,'RegionThreshold',1);
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
                        if norm(p0-tipping_pts{j}(1,:)) > norm(p1-tipping_pts{j}(1,:))
                            p_temp = p1; p1 = p0; p0 = p_temp;
                        end
                        p_dir = p1 - p0;p_dir = [p_dir(2),-p_dir(1)];
                        
         
                        
                        feature_bp{j} = mean(  sign((tipping_pts{j} - p0)*p_dir'));
                        
%                         if abs(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)) - max( tipping_pts{j}(:,1))+min(tipping_pts{j}(:,1))) < 0.01*alpha
%                             feature_bp{j} = 0;
%                         end
                        
%                         feature_bp{j} = tanh(10*feature_bp{j});
                        figure(11);
                        plot(shp_pos)
                        title( feature_bp{j} )
                        hold on,
                        plot([p0(1),p1(1)],[p0(2),p1(2)],'k','linewidth',2);
                        hold off;
                        end
                
                    
                case 'turbulence'
                    if length(tipping_times)>1000
                        tipping_times = tipping_times(end-1000:end);
                        tipping_pts{j} = tipping_pts{j}(end-1000:end,:);
                    end
                    if length(unique(tipping_times)) == length(tipping_times)
                        feature_bp{j} = 0;
                    else
                        time = mode(tipping_times);
                        idx = find( tipping_times == time);
                        feature_bp{j} = 1*( max( reshape(pdist( tipping_pts{j}(idx,:)' ),1,[]) )>10) ;
                    end
                    
                case 'area'
                    % this is only for Rossler
                    X = linspace(0,250,526); [X,Y] = meshgrid(X,X);
                    idx = find(u(:) > min(u(:))+0.9*(max(u(:)) - min(u(:))) );
                    shp = alphaShape( [X(idx),Y(idx)],1);
                    feature_bp{j} = sum( area(shp,1:numRegions(shp)))/250/250;
                    figure(11)
                    plot(shp)
                    title(feature_bp{j})
                case 'area2'
                    % this is only for Rossler
                    X = linspace(0,250,526); [X,Y] = meshgrid(X,X);
                    idx = find(u(:) > min(u(:))+0.9*(max(u(:)) - min(u(:))) );
                    shp = alphaShape( [X(idx),Y(idx)],1);
                    feature_bp{j} = area(shp,1:numRegions(shp))/250/250;
                    figure(11)
                    plot(shp)
                    title(feature_bp{j})
                case 'num'
                    % this is only for Rossler
                    X = linspace(0,250,526); [X,Y] = meshgrid(X,X);
                    idx = find(u(:) > min(u(:))+0.95*(max(u(:)) - min(u(:))) );
                    shp = alphaShape( [X(idx),Y(idx)],1);
                    feature_bp{j} = numRegions(shp);
                    figure(11)
                    plot(shp)
                    title(feature_bp{j})
                    
            end
                
                


            
        
    end
end

% fprintf('Ubp done\n')
if length(feature_bm) == 0
    for j = 1:Rep
        
        switch model
            case 'Barkley'
                anew = num2str(a-dn(1),'%1.6f');
                bnew = num2str(b-dn(2),'%1.6f');
                
                length(anew)
                length(bnew)
                
                fid = fopen('../DataGenerator/Spiral_Wave/Barkley/task.txt');
                C=textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                C{1}{1} = anew; C{1}{2} = bnew;
                writecell(C{1},'../DataGenerator/Spiral_Wave/Barkley/task.txt','QuoteStrings',0)
%                 if b > 0.05
%                     C{1}{11}(1:5) = num2str(100000);
%                 end
                cd ../DataGenerator/Spiral_Wave/Barkley
                system('rm tip.txt')
                system('touch tip.txt')
                system('./ezspiral')
                
                cd ../../../BifurcationTracing
                tipping_pts{j} = readtable('../DataGenerator/Spiral_Wave/Barkley/tip.txt');
                u = readtable('../DataGenerator/Spiral_Wave/Barkley/fc.txt');
                u = u.Var2;
                u = reshape(u(5:end),501,501);
                
            case 'Bar-Eiswirth'
                if a-dn(1) < 0
                    dn(1) = a-1E-14;
                end
                anew = num2str(1/(a-dn(1)),'%1.6f');
                bnew = num2str(b-dn(2),'%1.6f');
                
                length(anew)
                length(bnew)
                
               fid = fopen('../DataGenerator/Spiral_Wave/Bar-Eiswirth/task.txt');
                C=textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                C{1}{2} = bnew; C{1}{3} = anew;
                writecell(C{1},'../DataGenerator/Spiral_Wave/Bar-Eiswirth/task.txt','QuoteStrings',0)
                
                cd ../DataGenerator/Spiral_Wave/Bar-Eiswirth
                system('rm tip.txt')
                system('touch tip.txt')
                system('./ezspiral')
                
                cd ../../../BifurcationTracing/
                tipping_pts{j} = readtable('../DataGenerator/Spiral_Wave/Bar-Eiswirth/tip.txt');
                
                u = readtable('../DataGenerator/Spiral_Wave/Bar-Eiswirth/fc.txt');
                u = u.Var2;
                u = reshape(u(5:end),501,501);
                
            case 'Rossler'
                
                
                cnew = num2str(a-dn(1),'%1.6f');
                anew = num2str(b-dn(2),'%1.6f');
                
                fid = fopen('../DataGenerator/Spiral_Wave/Rossler/task.txt');
                C=textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                C{1}{3} = cnew; C{1}{1} = anew;
                writecell(C{1},'../DataGenerator/Spiral_Wave/Rossler/task.txt','QuoteStrings',0)

                cd ../DataGenerator/Spiral_Wave/Rossler

                system('./ezspiral')
                cd ../../../BifurcationTracing/

                u = readtable('../DataGenerator/Spiral_Wave/Rossler/fc.txt');
                u = u.Var2; 
                u = reshape(u(5:end),526,526);
%                 cd ../
                tipping_pts{j} = [];
        end
                if size(tipping_pts{j},1) == 0
                    tipping_pts{j} = [0,0];
                else
                    if strcmp(dist,'turbulence')
                        tipping_times = [tipping_pts{j}.Var1];
                    end
                        tipping_pts{j} = [tipping_pts{j}.Var2,tipping_pts{j}.Var3];
                    
                end
                switch model
                    case 'Barkley'
                       if size(tipping_pts{j},1) > 1000 || size(tipping_pts{j},1)==1
                            tipping_pts{j} = tipping_pts{j}(1000:end,:);

                        end
                        idx = find( (tipping_pts{j}(:,1)>98)+(tipping_pts{j}(:,1)<2) +...
                            + (tipping_pts{j}(:,2)>98) + (tipping_pts{j}(:,2)<2)  );
                        
                    case 'Bar-Eiswirth'
                        if strcmp(dist,'turbulence') == 0
                            if size(tipping_pts{j},1) > 2000 || size(tipping_pts{j},1)==1
                                tipping_pts{j} = tipping_pts{j}(2000:end,:);

                            end
                            idx = find( (tipping_pts{j}(:,1)>48)+(tipping_pts{j}(:,1)<2) +...
                                + (tipping_pts{j}(:,2)>48) + (tipping_pts{j}(:,2)<2)  );
                        else
                            idx = [];
                        
                        end
                        
                    case 'Rossler'
                        tipping_pts{j} = [];
                        idx = [];
                        
                end
                if length(idx)
                    tipping_pts{j} = tipping_pts{j}(1:idx(1),:);
                end
                if size(tipping_pts{j},1) && size(tipping_pts{j},2)>1
%                     idx = find( pdist2(tipping_pts{j} , tipping_pts{j}(1,:)) < 0.5 );
%                     idx = idx( find(idx>200));
%                     if length(idx)
%                         tipping_pts{j} = tipping_pts{j}(1:idx(end),:);
%                     end
                
                
                    figure(10);
                    scatter(tipping_pts{j}(:,1),tipping_pts{j}(:,2),[],1:length(tipping_pts{j}(:,1)))
                    hold on,
                    scatter(tipping_pts{j}([1,end],1),tipping_pts{j}([1,end],2),'k*')
                    hold off
                end
       
        
            switch dist
                case 'steady'
                    feature_bm{j} = min(size(tipping_pts{j},1)/1000,1);
                    feature_bm{j} = size(tipping_pts{j},1) < 1000;
                    feature_bm{j} = max(u(:))-min(u(:)) < 1E-2;
                case 'ring_fit'
                    if size(tipping_pts{j},1) < 100
                        feature_bm{j} = 0;
                    else
%                         alpha0 = max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2));
%                         alpha0 = max( max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)),alpha0);
%                         alpha = alpha0/2.5;
%                         shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha,'RegionThreshold',1);
%                         while numRegions(shp_pos) < 1
%                             alpha = alpha*0.9;
%                             shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha);%,alpha,'RegionThreshold',1);
%                             
%                             
%                         end
%                         shp_area = area(shp_pos,1);
%                         shp_peri = perimeter(shp_pos,1)/2/pi;
% %                         shp_rdns = shp_area/(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)));
% %                         shp_rdns = shp_rdns/(max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)));%/pi/shp_peri/shp_peri;
% %                         shp_area = shp_area/(   sum(perimeter(shp_pos,1:numRegions(shp_pos)))^2);
%                         shp_rdns = shp_area/pi/shp_peri/shp_peri;

                        alpha0 = max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2));
                        alpha0 = max( max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)),alpha0);
                        alpha = alpha0/2.5;
%                         shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha,'RegionThreshold',1,'HoleThreshold',1);
%                         while numRegions(shp_pos) < 1
%                             alpha = alpha*0.9;
%                             shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha);%,alpha,'RegionThreshold',1);
%                             
%                             
%                         end
%                         shp_area = area(shp_pos,1);
%                         shp_peri = perimeter(shp_pos,1)/2/pi;
% %                         shp_rdns = shp_area/(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)));
% %                         shp_rdns = shp_rdns/(max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)));
%                         shp_rdns = shp_area/pi/shp_peri/shp_peri;
                        
                        shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),alpha);%,alpha,'RegionThreshold',1);
                        
                        shp_area = sum(area(shp_pos,1:numRegions(shp_pos)));
                        shp_rdns = shp_area/(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)));
                        shp_rdns = shp_rdns/(max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)));
%                         
                        feature_bm{j} = tanh(4*shp_rdns);
                        
                        figure(11);
                        plot(shp_pos)
                        title( feature_bm{j} )
                    end
                        
                    
                case 'line_fit'
                    if size(tipping_pts{j},1) < 100
                          feature_bm{j} = 0;
                     else
                        out = eigs( ( tipping_pts{j}- mean(tipping_pts{j}))'*(tipping_pts{j} - mean(tipping_pts{j})));
                        l = out(2)/out(1);
                        l
                        while l > 0.6
                            tipping_pts{j} = tipping_pts{j}(1:round(size(tipping_pts{j},1)*0.75),:);
                            out = eigs( ( tipping_pts{j}- mean(tipping_pts{j}))'*(tipping_pts{j} - mean(tipping_pts{j})));
                            l = out(2)/out(1);
                        end
                        
                        alpha = max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2));
                        alpha = max( max( tipping_pts{j}(:,1))-min(tipping_pts{j}(:,1)),alpha);
                        
                        
                        
                        
                        shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2),10*alpha,'RegionThreshold',1);
                        if numRegions(shp_pos) < 1
                            
                            shp_pos = alphaShape(tipping_pts{j}(:,1),tipping_pts{j}(:,2));%,alpha,'RegionThreshold',1);
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
                        if norm(p0-tipping_pts{j}(1,:)) > norm(p1-tipping_pts{j}(1,:))
                            p_temp = p1; p1 = p0; p0 = p_temp;
                        end
                        p_dir = p1 - p0; p_dir = [p_dir(2),-p_dir(1)];
                        feature_bm{j} = mean(  sign((tipping_pts{j} - p0)*p_dir'));
                        
                        
%                         if abs(max( tipping_pts{j}(:,2))-min(tipping_pts{j}(:,2)) - max( tipping_pts{j}(:,1))+min(tipping_pts{j}(:,1))) < 0.01*alpha
%                             feature_bm{j} = 0;
%                         end
                        
                        figure(11);
                        plot(shp_pos)
                        title( feature_bm{j} )
                        hold on,
                        plot([p0(1),p1(1)],[p0(2),p1(2)],'k','linewidth',2);
                        hold off;
                    end
                    
               
                case 'turbulence'
                    if length(tipping_times)>1000
                        tipping_times = tipping_times(end-1000:end);
                        tipping_pts{j} = tipping_pts{j}(end-1000:end,:);
                    end
                    if length(unique(tipping_times)) == length(tipping_times)
                        feature_bm{j} = 0;
                    else
                        time = mode(tipping_times);
                        idx = find( tipping_times == time);
                        feature_bm{j} = 1*( max( reshape(pdist( tipping_pts{j}(idx,:)' ),1,[]) )>10) ;
                    end
                    
                case 'area'
                    % this is only for Rossler
                    X = linspace(0,250,526); [X,Y] = meshgrid(X,X);
                    idx = find(u(:) > min(u(:))+0.9*(max(u(:)) - min(u(:))) );
                    shp = alphaShape( [X(idx),Y(idx)],1);
                    feature_bm{j} = sum( area(shp,1:numRegions(shp)))/250/250;
                    figure(11)
                    plot(shp)
                    title(feature_bm{j})
                case 'area2'
                    % this is only for Rossler
                    X = linspace(0,250,526); [X,Y] = meshgrid(X,X);
                    idx = find(u(:) > min(u(:))+0.9*(max(u(:)) - min(u(:))) );
                    shp = alphaShape( [X(idx),Y(idx)],1);
                    feature_bm{j} = area(shp,1:numRegions(shp))/250/250;
                    figure(11)
                    plot(shp)
                    title(feature_bm{j})
                case 'num'
                    % this is only for Rossler
                    X = linspace(0,250,526); [X,Y] = meshgrid(X,X);
                    idx = find(u(:) > min(u(:))+0.95*(max(u(:)) - min(u(:))) );
                    shp = alphaShape( [X(idx),Y(idx)],1);
                    feature_bm{j} = numRegions(shp);
                    figure(11)
                    plot(shp)
                    title(feature_bm{j})
            
            end
                
    end
end

if strcmp(dist,'area2')
    out = ws_distance(feature_bp{1}, feature_bm{1}, 2);
else
    out = abs( feature_bp{1} - feature_bm{1} );
end
% fprintf('Ubm done\n')

