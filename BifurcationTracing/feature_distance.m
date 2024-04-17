function [out,H] = feature_distance(feature_bp_all, feature_bm_all, dist)


out = 0; H = 0;
switch dist
    
        
   
    case 'num'

        out = out + ws_distance( feature_bm_all,feature_bp_all, 2 );

    case 'roundness-bag'
                
        out = out + ws_distance( feature_bm_all,feature_bp_all, 2 );
               
        
        
    case 'roundness'
        Rep = length(feature_bp_all);
        dist_mat = zeros( Rep, Rep);

        for i = 1:Rep
            for j = 1:Rep
%                             C = pdist2(feature_bm_all{i}', feature_bp_all{j}');
%                             dist_mat(i,j) = mexEMD( ones(length(feature_bm_all{i}),1)/length(feature_bm_all{i}) , ones(length(feature_bp_all{j}),1)/length(feature_bp_all{j}), C);
                dist_mat(i,j) = ws_distance( feature_bm_all{i},feature_bp_all{j}, 2 );


            end
        end
%         cd mexEMD
        out = mexEMD( ones(Rep,1)/Rep, ones(Rep,1)/Rep, dist_mat);
%         cd ..

    otherwise
        out = abs( feature_bp_all{1} - feature_bm_all{1} );



        
end