function [R,b,x1,x2,theta,p1,p2,py] = find_rotation_translation_alpha(U1, U2, X, Y, tip1, tip2, r, h, init)


% This function returns the rotation matrix R and translation vector b from
% the spiral pattern in U1 to U2. 
% Relation: x2 = R*x1 + b, or x1 = R'*(x2 - b)

if nargin < 9
    init = [0,0,0]';
end




idx_cut = find(  (X(:) - tip2(1)).^2 + (Y(:) - tip2(2)).^2 > r^2);
U2( idx_cut ) = min(U2(:)) + 0.1*( max(U2(:)) - min(U2(:)));


idx = find(  U2 > 0.5*( max(U2(:)) - min( U2(:))) + min( U2(:))     );
x2 = [X(idx),Y(idx)];
U2 = U2(idx);

shp2 = alphaShape(x2(:,1),x2(:,2),h);

if numRegions(shp2) > 1
    areas = area(shp2, 1:numRegions(shp2));
    idx = find(areas == max(areas));
else
    idx = 1;
end

[~,P2] = boundaryFacets(shp2,idx);
p2 = polyshape(P2(:,1)',P2(:,2)');




% idx_cut = find(  (X(:) - tip1(1)).^2 + (Y(:) - tip1(2)).^2 > r^2);
% U1( idx_cut ) = min(U1(:)) + 0.1*( max(U1(:)) - min(U1(:)) );


idx = find(  U1 > 0.5*( max(U1(:)) - min( U1(:))) + min( U1(:))     );
x1 = [X(idx),Y(idx)];
U1 = U1(idx);




shp1 = alphaShape(x1(:,1),x1(:,2),h);

if numRegions(shp1) > 1
    areas = area(shp1, 1:numRegions(shp1));
    idx = find(areas == max(areas));
else
    idx = 1;
end

[~,P1] = boundaryFacets(shp1,idx);
p1 = simplify(polyshape(P1(:,1)',P1(:,2)'));




Translate = @(theta, b,x) [cos(theta), -sin(theta); sin(theta), cos(theta)]*(x' + b);

obj = @(theta, b) 0*norm( Translate(theta,b,tip1') - tip2 ) - Area_alpha(theta,b,P1,p2);  %1/size_c*norm( Translate(theta,b,x1) - x2');
options = optimset('MaxFunEvals',1E6);
out = fminsearch( @(x) obj(x(1),x(2:3)), init,options);

theta = out(1); 
b = out(2:3);
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

[~,py] = Area_alpha(theta,b,P1,p2);

function [a,py] = Area_alpha(theta, b, P1, p2)

% y = Translate(theta,b,x1)';
% shpy = alphaShape(y(:,1),y(:,2),h);
% [~,Py] = boundaryFacets(shpy,1);

y = Translate(theta,b,P1)';

% idx_cut = find(  (y(:,1) - tip2(1)).^2 + (y(:,2) - tip2(2)).^2 > r^2);
% y = y(idx_cut, :);

py = simplify(polyshape(y(:,1)',y(:,2)'));


% py = polyshape(Py(:,1)',Py(:,2)');
a = area(simplify(intersect(py,p2)));
end

end
    
    