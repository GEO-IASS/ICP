%% load data and randomly generate transformation
folder = './drill';

%load model 
load(strcat(folder, '/model.mat'))
figure(1)
P_m = Mdata;

%randomly rotate the cloud
eul = rand(1, 3)*2*pi;
R = eul2rotm(eul);
%R = eul2rotm([0, 0, 0.1]);
P_n = R*P_m;

[nd, np] = size(P_m);

%randomly translate the cloud
T = rand(3, 1)*10 - 5;
P_n = bsxfun(@plus, P_n, T);
pcshow(P_n')

%% ICP initialization
%initialize hypothesis on R and T
R_est = eye(3);

%given the hypothesis, the proposed transformed point cloud
P_n_est = P_m;
hold on
pcshow(P_n_est');

%construct KD tree on the measured points
Md_tree = KDTreeSearcher(P_n');
error = norm(P_n_est - P_n, 'fro')

%% Iterative Closest Point
%knnsearch
while error > 1e-5*np
    
    knnind = knnsearch(Md_tree, P_n_est');
    
    %first calculate distance between registered points
    dist = sum((P_n_est' - P_n(:, knnind)').^2, 2);
    
    %then find the median, reject points that are too far away(below the median)
    %indeff = dist > median(dist);
    
    %here trimmed points
    %P_n_trim = P_n(:, knnind(indeff))';
    P_n_trim = P_n';
    %P_m_trim = P_m(:, indeff);
    P_m_trim = P_m;
    
    %do SVD to estimate the transformation
    %   first remove the translation
    C_n_trim = mean(P_n_trim, 1);
    C_m_trim = mean(P_m_trim, 2);
    T_est = C_n_trim' - C_m_trim;
    
    
    Q = bsxfun(@minus, P_m_trim ,C_m_trim);
    Q2 = bsxfun(@minus, P_n_trim, C_n_trim);
    %   calculate H matrix
    H = Q*Q2;
    %   call SVD and follow the formula in the paper fool
    [U, S, V] = svd(H);
    X = V*U';
    
    if det(X) > 0
        R_est = X;
        P_n_est = bsxfun(@plus, R_est*P_m, T_est);
        error = norm(P_n_est - P_n, 'fro')
    else
        disp('nooooo')
    end
end

figure(2)
pcshow(P_n_est')
hold on
pcshow(P_n')
