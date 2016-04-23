    figure(2)
    pcshow(P_model)
    hold on
    pcshow(P_rough, C_)
    hold off
    pause(0.1)
    %{
    P_proposed = bsxfun(@plus, P_rough*R_fine', T_fine');
    figure(2)
    pcshow(P_model)
    hold on
    pcshow(P_proposed)
    hold off
    %}