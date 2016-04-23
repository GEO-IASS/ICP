function [P_, C_] = remove_floor(P, C, iterations, episilon)
    %REMOVE_FLOOR removes floor using ransac to remove dominant plane
    %input P is point cloud, iteration is the number of total iterations
    %for ransac. Output argument ind is the indicies of the points which
    %have the floor removed
    
    best_outliers = inf;  
    
    for i = 1:iterations
        %pick 3 random points
        PointsPicked = datasample(P, 3, 1, 'Replace', false);
        planevec = PointsPicked\[1;1;1];
        
        DisFromPlane = (abs(P*planevec - 1))/norm(planevec);
        outliers = (DisFromPlane > episilon);
        num_outliers = sum(outliers);
        
        if num_outliers < best_outliers
            ind = outliers;
            best_outliers = num_outliers;
        end
    end
    %select ponts
    P_ = P(ind, :);
    C_ = C(ind, :);
    
    %find the center of gravity and clean up anything that is too far out
    for i = 1:2
        P_center = mean(P_, 1);
        d_c = sum((bsxfun(@minus, P_, P_center).^2), 2);
        ind = (d_c) < 3*median(d_c);
        
        P_ = P_(ind, :);
        C_ = C_(ind, :);
    end
end

