function P_sampled = NS( P, np )
    %P is point cloud, np: number of points desired sample
    %resolution the sphere is split up into 100 bins
    
    PC_obj = pointCloud(P);
    normals = pcnormals(PC_obj);
    
    %put into bins
    % z = cos(phi)
    nu = (normals(:, 3) + 1)/2;
    theta = atan2(normals(:, 2), normals(:, 1));
    theta(theta<0) = 2*pi + theta(theta<0);
    
    nu_bin = ceil(nu/0.1) - 1;
    theta_bin = ceil(theta/2/pi*10) - 1;
    
    bin = (nu_bin*10 + theta_bin);
    
    bin_avai = unique(bin);
    np_per_bin = floor(np/length(bin_avai));
    %sample uniformaly in normal (spac
    %for each bin, sample np points, if the bin is empty, skip bin
    P_sampled = zeros(np, 3);
    count = 1;
    for i = 1:length(bin_avai)
        b = bin_avai(i);
        start = 1 + (count-1)*np_per_bin;
        finish = count*np_per_bin;
        P_bin = P(bin == b, :);
        P_sampled(start:finish, :)= datasample(P_bin, np_per_bin, 1);
        count = count + 1;
    end
    P_sampled(finish+1:end, :) = [];
end
    
