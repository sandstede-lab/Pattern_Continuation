function [out,H] = feature_distance(feature_bp_all, feature_bm_all, dist, mode)



switch mode
    case 'ks'
        switch dist
            case 'num'
                feature_bp_all = cat(1,feature_bp_all{:});
                feature_bm_all = cat(1,feature_bm_all{:});
                [H,out] = ttest2(feature_bp_all, feature_bm_all);
            case 'roundness'
                [~,~,out] = kstest2(feature_bp_all, feature_bm_all);
        end
        
        
    case 'dist'
        
        out = 0; H = 0;
        switch dist
            case 'num'
                
                
%                 for i = 1:Rep
%                     for j = 1:Rep
%                         out = out + abs( feature_bp_all{i} - feature_bm_all{j} );
% 
%                     end
%                 end
% 
% 
%         %         out
%                 out = out/Rep/Rep;
                
%                     feature_bp_all = cat(1,feature_bp_all{:});
%                     feature_bm_all = cat(1,feature_bm_all{:});

                    out = out + abs(sum(feature_bp_all) - sum(feature_bm_all))/length(feature_bp_all);
                

            otherwise
                
                    out = out + ws_distance( feature_bm_all,feature_bp_all, 2 );
               
        end
        
    case 'metric'
        out = 0; H = 0;
        switch dist
            case 'num'
                
                
%                     Rep = length(feature_bp_all);
%                     dist_mat = zeros( Rep, Rep);
%                     for i = 1:Rep
%                         for j = 1:Rep
%                             dist_mat(i,j) = abs(feature_bm_all{i}-feature_bp_all{j});
%                             
%                         end
%                     end
%                     cd mexEMD
%                     out = mexEMD( ones(Rep,1)/Rep, ones(Rep,1)/Rep, dist_mat);
%                     cd ..
                    out = out + ws_distance( feature_bm_all,feature_bp_all, 2 );
               
            otherwise
                    Rep = length(feature_bp_all);
                    dist_mat = zeros( Rep, Rep);
                    
                    for i = 1:Rep
                        for j = 1:Rep
%                             C = pdist2(feature_bm_all{i}', feature_bp_all{j}');
%                             dist_mat(i,j) = mexEMD( ones(length(feature_bm_all{i}),1)/length(feature_bm_all{i}) , ones(length(feature_bp_all{j}),1)/length(feature_bp_all{j}), C);
                            dist_mat(i,j) = ws_distance( feature_bm_all{i},feature_bp_all{j}, 2 );
                            
                            
                        end
                    end
                    cd mexEMD
                    out = mexEMD( ones(Rep,1)/Rep, ones(Rep,1)/Rep, dist_mat);
                    cd ..
                    
%                     feature_bm_all_list = zeros(1,Rep); feature_bp_all_list = zeros(1,Rep);
%                     for i = 1:Rep
%                         feature_bm_all_list(i) = mean(feature_bm_all{i});
%                         feature_bp_all_list(i) = mean(feature_bp_all{i});
%                     end
%                     out = out + ws_distance( feature_bm_all_list,feature_bp_all_list, 2 );
                 
                    
               
        end
        
end