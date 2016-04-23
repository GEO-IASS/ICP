close all
T_hist_wm = zeros(nframes, 3);
q_hist_wm = zeros(nframes, 4);

for i = 1:nframes
    q_c_to_w = pose(i, 3:6);
    q_c_to_m = q_hist(i, :);
    
    T_c_to_w = pose(i, 7:9);
    T_c_to_m = T_hist(i, :);
    
    R_c_to_w = quat2rotm(q_c_to_w);
    R_c_to_m = quat2rotm(q_c_to_m);
    R_w_to_c = R_c_to_w';
    
    R_w_to_m = R_c_to_m*R_w_to_c;
    T_w_to_m = T_c_to_m - (R_c_to_m*R_w_to_c*T_c_to_w')'; 
    
    T_hist_wm(i, :) = T_w_to_m;
    q_hist_wm(i, :) = rotm2quat(R_w_to_m);
    
end
   figure(1)
   plot(q_hist(:, 1))
   hold on
   plot(q_hist_wm(:, 1))
   plot(pose(:,3))
   legend('c2m', 'w2m', 'c2w')
   title('q1')
   
   figure(2)
   plot(q_hist(:, 2))
   hold on
   plot(q_hist_wm(:, 2))
   plot(pose(:,4))
   legend('c2m', 'w2m', 'c2w')
   title('q2')
   
   figure(3)
   plot(q_hist(:, 3))
   hold on
   plot(q_hist_wm(:, 3))
   plot(pose(:,5))
   legend('c2m', 'w2m', 'c2w')
   title('q3')
   
   figure(4)
   plot(q_hist(:, 4))
   hold on
   plot(q_hist_wm(:, 4))
   plot(pose(:,6))
   legend('c2m', 'w2m', 'c2w')
   title('q4')
   
   figure(5)
   plot(T_hist(:,1))
   hold on
   plot(T_hist_wm(:, 1))
   plot(pose(:, 7))
   legend('c2m', 'w2m', 'c2w')
   title('x')
   
   figure(6)
   plot(T_hist(:,2))
   hold on
   plot(T_hist_wm(:, 2))
   plot(pose(:, 8))
   legend('c2m', 'w2m', 'c2w')
   title('y')
   
   figure(7)
   plot(T_hist(:,3))
   hold on
   plot(T_hist_wm(:, 3))
   plot(pose(:, 9))
   legend('c2m', 'w2m', 'c2w')
   title('z')
   
   