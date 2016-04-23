close all;
clear;

%% Initialization (load model..things that I only do once)

%  sensor calibration parameters
addpath('./params');
load calib_vicon_xtion.mat            % s_cam3, R_cv3, t_cv3 camera-vicon calibration parameters

folder = './liq_container';

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
axang_episilon = 0.002;
mean_error_thresh = 0.2;
num_ipose = 10;

%things that depend on the parameters
grid729 = zeros(3, 27);
pos = [0, 50, -50];
ori = [0, 90, -90];

count = 1;
for i = 1:3
    for j = 1:3
        for k = 1:3
            grid729(1:3, count) = [pos(i), pos(j), pos(k)];
            count = count + 1;
        end
    end
end
         

%initialization
drill_tree = KDTreeSearcher(P_model); %kd tree for the model

%% For all the frames

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
    %R_CW = quat2rotm(qcam);
    %P_ = bsxfun(@plus, R_CW * P_', xcam')';
    
    %colors for visualization
    C_ = [r g b]./255;
    
         figure(1),
         pcshow(P_,[r g b]/255)
    
    %find and remove floor
    [P_, C_] = remove_floor(P_, C_, 1000, episilonRANSAC);
    
    figure(5)
    pcshow(P_, C_)
    
    pause(0.1)
    
    
    %     figure(2)
    %     pcshow(P_, C_)
    %
    %     pause(1)
    
    %%%%%% Generate proper initial pose and run ICP (possibly multiple times) to estimate the pose
    % you may use KDTreeSearcher and knnsearch for correspondence search.
    
    error_best = inf;
    
    
    [npmeas, ~] = size(P_);
    ctr_ = mean(P_, 1);
    Q_ = bsxfun(@minus, P_, ctr_)';
    ctr_ = ctr_';
    
    if k == 1
        %rough estimate just by shifting
        R_rough = eye(3);
        %T_rough = [0; 0; 0];
        T_rough = mean(P_model, 1)' - ctr_;
    else
        R_rough = R_est;
        T_rough = T_est;
    end
    
    P_rough = bsxfun(@plus, P_*R_rough', T_rough');
    showrough
    
    
    
    %normal space sampling from P_rough to match
    [P_rough, C_] = NS_showy(P_rough, C_, (ceil(npmeas/2)));
    ctr_rough = mean(P_rough, 1);
    Q_rough = bsxfun(@minus, P_rough, ctr_rough)';
    ctr_rough = ctr_rough';
    [npmeas, ~] = size(P_rough);
    showrough2
    
    for i = 1:27
        
        R_fine = eye(3);
        T_fine = grid729(1:3, i);
        
        %do one iteration of ICP given the random jitter
        R_prev = R_fine;
        T_prev = T_fine;
        axang_diff = 10;
        T_diff = 10;
        
        while axang_diff > axang_episilon && T_diff > 0.1
            
            [R_fine, T_fine, detX, error] = ICP_showy(P_rough, P_model, Q_rough, ctr_rough, drill_tree, R_fine, T_fine, C_);
            
            if detX < 0
                fprintf('determinant of R is -1, skip this iteration')
                break;
            end
            
            %match and calculate error again
            axang = rotm2axang(R_fine*R_prev');
            
            %print error
            if isempty(notetext)
                errortext = sprintf('%f', error);
                notetext = sprintf('search pose %d, error is %s', i, errortext);
                l_e_prev = length(errortext);
                fprintf(notetext)
            else
                errortext = sprintf('%f, %d', error);
                backspaces = repmat(sprintf('\b'), 1, l_e_prev);
                fprintf(backspaces)
                fprintf(errortext)
                
                l_e_prev = length(errortext);
            end
            
            axang_diff = axang(4);
            T_diff = norm(T_fine - T_prev);
            R_prev = R_fine;
            T_prev = T_fine;
        end
        
        fprintf('\n')
        notetext = [];
        
        if error < error_best
            error_best = error;
            meanerror = sqrt(error/npmeas);
            R_best = R_fine;
            T_best = T_fine;
            P_proposed = bsxfun(@plus, P_rough*R_fine', T_fine');
            if meanerror < mean_error_thresh
                R_est = (R_best*R_rough);
                T_est = R_best*T_rough + T_best;
                % log results
                R_record = R_est;
                T_record = T_est;
                q_hist(k, :) = rotm2quat(R_record);
                T_hist(k, :) = T_record';
                break;
            end
        end
        R_est = (R_best*R_rough);
        T_est = R_best*T_rough + T_best;
        % log results
        R_record = R_est;
        T_record = T_est;
        q_hist(k, :) = rotm2quat(R_record);
        T_hist(k, :) = T_record';
        
    end
    
end




