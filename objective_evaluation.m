function out = objective_evaluation(feature1, feature2)

% This function outputs distance between feature collections



out = 0; 

type = checkType(feature1);

switch type
    
        
   
    case 'scalar'

        out = abs( feature1 - feature2 );

    case 'vector'
                
        out = ws_distance( feature1, feature2, 2 );
               
        
        
    case 'cell'
        path(path,'mexEMD/');
        Rep = length(feature1);
        dist_mat = zeros( Rep, Rep);

        for i = 1:Rep
            for j = 1:Rep
%                             C = pdist2(feature_bm_all{i}', feature_bp_all{j}');
%                             dist_mat(i,j) = mexEMD( ones(length(feature_bm_all{i}),1)/length(feature_bm_all{i}) , ones(length(feature_bp_all{j}),1)/length(feature_bp_all{j}), C);
                dist_mat(i,j) = ws_distance( feature1{i},feature2{j}, 2 );


            end
        end
%         cd mexEMD
        out = mexEMD( ones(Rep,1)/Rep, ones(Rep,1)/Rep, dist_mat);
%         cd ..

    
end

function type = checkType(input)


if iscell(input)
    type = 'cell';
elseif length(input)==1
    type = 'scalar';
elseif isvector(input)
    type = 'vector';



end