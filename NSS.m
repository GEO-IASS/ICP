close all;
clear;

%% Initialization (load model..things that I only do once)

%  sensor calibration parameters
addpath('./params');
load calib_vicon_xtion.mat            % s_cam3, R_cv3, t_cv3 camera-vicon calibration parameters

folder = './drill';

% Load model
load(strcat(folder,'/model.mat'));
P_model = Mdata';

% Load cam pose
load(strcat(folder,'/pose.mat'))
[nframes, ~] = size(pose);
notetext = [];

q_hist = zeros(nframes, 4);
T_hist = zeros(nframes, 3);

%parameters
episilonRANSAC = 35;
axang_episilon = 0.001;
mean_error_thresh = 0.03;

%% Do many times

for k = 1:nframes
    fprintf('frame %d\n', k)
    
    frame = pose(k,1);
    t = pose(k,2);
    qcam = pose(k,3:6); % camera model orientation w.r.t. reference coordinate frame
    xcam = pose(k,7:9); % camera model position w.r.t. reference coordinate frame
    
    % Load images
    depth = imread(strcat(folder,'/depth/',sprintf('%06d.png',k)));
    rgb = imread(strcat(folder,'/rgb/',sprintf('%06d.png',k)));
    
    % point cloud
    [pcx, pcy, pcz, r, g, b, D_, X, Y,validInd] = depthToCloud_full_RGB(depth, rgb, './params/calib_xtion.mat');
    % NOTE:
    % - pcz equals to D_(validInd)
    % - pcx equals to X(validInd)
    % - pcy equals to Y(validInd)
    
    % Adjust actual camera to camera model.
    % (Please don't worry too much about difference between "actualy camera" and
    % "camera model". Just do not remove the line below.)
    P_ = s_cam3*[pcx pcy pcz]*[0 0 1;-1 0 0;0 -1 0]'*R_cv3 + repmat(t_cv3,length(pcx),1);
    %convert measured points in to world frame
    R_CW = quat2rotm(qcam);
    P_ = bsxfun(@plus, R_CW * P_', xcam')';
    
%     figure(1),
%     pcshow(P_,[r g b]/255)
    
    
    
    %%%%%% Find and remove the floor
    
    numpts = length(pcx);
    best_inliers = 0;
    lnliers_ind = 0;
    
    drill_radius = 500;
    
    
    %do ransac 500 times
    for i = 1:500
        %pick 3 random points
        picks = randsample(numpts, 3);
        PointsPicked = P_(picks, :);
        planevec = PointsPicked\[1;1;1];
        
        DisFromPlane = (abs(P_*planevec - 1))/norm(planevec);
        inliers = (DisFromPlane < episilonRANSAC);
        num_inliers = sum(inliers);
        
        if num_inliers > best_inliers
            inliers_ind = inliers;
            bestvec = planevec;
            best_inliers = num_inliers;
        end
    end
    
    ind = ~inliers_ind;
    P_meas = P_(ind, :);
    
    colors = [r(ind) g(ind) b(ind)]./255;
    
    %find the center of gravity and clean up anything that is too far out
    P_center = mean(P_meas, 1);
    d_c = sum((bsxfun(@minus, P_meas, P_center).^2), 2);
    ind = (d_c - mean(d_c)) < 2*std(d_c);
    
    % pre calculate this for ICP
    P_meas = P_meas(ind, :);
    C_meas = mean(P_meas, 1);
    Q_meas = bsxfun(@minus, P_meas, C_meas)';
    %
    
    colors = colors(ind, :);
%     %figure(2)
%     pcshow(P_meas, colors)
%     
%     pause(1)
    
    %%%%%% Generate proper initial pose and run ICP (possibly multiple times) to estimate the pose
    % you may use KDTreeSearcher and knnsearch for correspondence search.
    
    if k == 1
        
        % space to search for EM
        num_ipose = 10;
        num_zoff = 10;
        yaws = linspace(0, 2*pi, num_ipose);
        zoffs = linspace(-10, 10, num_zoff);
        
        %initialize parameters that stores the best results
        error_best = inf;
        R_best = eye(3);
        T_best = zeros(3, 1);
        term_early = 0; % flag that signify whether i can terminate early
        
        [npmeas, ~] = size(P_meas); %number of measured points (after floor removal)
        
        %construct KD tree of the model
        disp('Constructing KD tree from model point cloud...\n')
        drill_tree = KDTreeSearcher(P_model);
        
        %calculate normal for each point in the model
        disp('Constructing matlab ptcloud model...')
        PO_model = pointCloud(P_model);
        PO_meas = pointCloud(P_meas);
        disp('Calculating normal vectors for each point in the model...\n')
        N_model = pcnormals(PO_model);
        N_meas = pcnormals(PO_meas);
        
        %bin normals
        disp('Binning Normal...10 bins for both polar and azimuth')
        
        %model
        phi_model = acos(N_model(:, 3));
        theta_model = atan2(N_model(:,2), N_model(:,1));
        bin_model = ceil(10*((phi_model)./(pi/10) - 1) + (theta_model+pi)./(2*pi/10));
        %bin_model_cell = cell(100,1);
        
        %measured
        phi_meas = acos(N_meas(:, 3));
        theta_meas = atan2(N_meas(:,2), N_meas(:,1));
        bin_meas = ceil(10*((phi_meas)./(pi/10) - 1) + (theta_meas+pi)./(2*pi/10));
        %bin_meas_cell = cell(100,1);
        
        %randomly sample from each bin now let's do 2000 samples
        num_sample = 2000;
        sample_model = zeros(num_sample, 3);
        sample_meas = zeros(num_sample, 3);
        
        
        for b = 1:100
            sample_model = datasample(P_model(bin_model == b, :), num_sample/100, 1);
            sample_meas = datasample(P_model(bin_meas == b, :), num_sample/100, 1);
            %bin_model_cell{b} = P_model(bin_model == b, :);
            %bin_meas_cell{b} = P_meas(bin_meas == b, :);
        end
        
        
        %Do normal space sampling to get a good starting point
        
        
        %get initial pose by just shifting the measured points
        %{
        eul = [0, 0, 0];
        R_rough = eul2rotm(eul);
        P_rough = (R_rough*Q_meas)';
        
        %find nearest neighbors
        ind = knnsearch(drill_tree, P_rough);
        P_match = P_model(ind, :);
        %calculate error
        error = norm(P_match - P_rough, 'fro');
        error_prev = 0;
        
        C_model = mean(P_match, 1);
        T_rough = (C_model - C_meas)';
        P_rough = bsxfun(@plus, P_meas, T_rough')';
        
        C_rough = mean(P_rough, 2);
        Q_rough = bsxfun(@minus, P_rough, C_rough);
        %}
    else
        %use the results from last iteration
        error_best = inf;
        R_best = eye(3);
        T_best = zeros(3, 1);
        R_rough = R_est';
        T_rough = -R_rough*T_est;
        term_early = 0;
        [npmeas, ~] = size(P_meas);
        P_rough = bsxfun(@plus, R_est'*P_meas', T_rough);%initialize with last steps stuff
        C_rough = mean(P_rough, 2);
        Q_rough = bsxfun(@minus, P_rough, C_rough);
    end
    
        
        

    %now do fine estimation with different initializations
    for i = 1:num_ipose
        for j = 1:num_zoff
            eul = [yaws(i), 0, 0];
            R_fine = eul2rotm(eul);
            T_fine = rand(3,1)*20;
            %T_fine = [0 0 0]';
            C_fine = C_rough + T_fine;
            P_proposed = bsxfun(@plus, (R_fine*P_rough), T_fine)';
            %find nearest neighbors
            ind = knnsearch(drill_tree, P_proposed);
            P_match = P_model(ind, :);
            %calculate error
            error = norm(P_match - P_proposed, 'fro');
            error_prev = 0;
            R_prev = R_fine;
            axang_diff = 10;
            
            while axang_diff > axang_episilon
                %ICP here
                error_prev = error;
                
                %estimate translation
                C_model = mean(P_match, 1);
                
                %estimate rotation
                Q_model = bsxfun(@minus, P_match, C_model);
                H = Q_rough*Q_model;
                [U, S, V] = svd(H);
                X = V*U';
                
                %update proposal
                R_fine = X;
                T_fine = C_model' - R_fine*C_rough;
                P_proposed = bsxfun(@plus, R_fine*P_rough, T_fine)';
                
                %match and calculate error again
                error = norm(P_match - P_proposed, 'fro');
                axang = rotm2axang(R_fine*R_prev');
                if isempty(notetext)
                    errortext = sprintf('%f', error);
                    notetext = sprintf('search pose %d, position %d, error is %s', i, j, errortext);
                    l_e_prev = length(errortext);
                    fprintf(notetext)
                else
                    errortext = sprintf('%f', error);
                    backspaces = repmat(sprintf('\b'), 1, l_e_prev);
                    fprintf(backspaces)
                    fprintf(errortext)
                    l_e_prev = length(errortext);
                end
                axang_diff = axang(4);
                R_prev = R_fine;
                ind = knnsearch(drill_tree, P_proposed);
                P_match = P_model(ind, :);
            end
            
            fprintf('\n')
            notetext = [];
            
            if error < error_best
                error_best = error;
                meanerror = error/npmeas;
                R_best = R_fine;
                T_best = T_fine;
                if meanerror < mean_error_thresh
                    term_early = 1;
                    break;
                end
            end
        end
        if term_early
            break;
        end
    end
    
    %%%%%% Show the estimated pose in the reference frame
    R_est = (R_best*R_rough)';
    T_est = -(R_rough'*T_rough + R_est*T_best);
    % log results
    q_hist(k, :) = rotm2quat(R_est);
    T_hist(k, :) = T_est';
    
    %P_proposed = bsxfun(@plus, R_fine_best*P_rough, T_fine_best)';
    
    P_model_t = bsxfun(@plus, R_est*P_model', T_est)';
    
%     figure(3)
%     pcshow(P_model)
%     hold on
%     pcshow(P_model_t)
%     pcshow(P_meas, colors)
%     hold off   
end
    




