function pointLiesInside = icatb_check_point_inside(currentPoint, axesLocation)
% Use cross product to determine if the point lies inside
 
vertices = zeros(5, 2);
vertices(1, :) = [axesLocation(1), axesLocation(2)];
vertices(2, :) = [axesLocation(1) + axesLocation(3), axesLocation(2)];
vertices(3, :) = [axesLocation(1) + axesLocation(3), axesLocation(2) + axesLocation(4)];
vertices(4, :) = [axesLocation(1), axesLocation(2) + axesLocation(4)];
vertices(5, :) = vertices(1, :);
 
cross_product = zeros(1, size(vertices, 1) - 1);
for kk = 1:size(vertices, 1) - 1
    vectorA = vertices(kk + 1, :) - vertices(kk, :);
    vectorB = currentPoint - vertices(kk, :);
    cross_product = vectorA(1, 1)*vectorB(1, 2) - vectorA(1, 2)*vectorB(1, 1);
    if kk == 1
        pointLiesInside =  (cross_product >= 0);
    else
        pointLiesInside = pointLiesInside & (cross_product >= 0);
    end
end
