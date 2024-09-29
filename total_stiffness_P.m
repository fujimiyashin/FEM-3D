function res = total_stiffness_P(unit_stiffness, node_sum, node_id_lis)
    Ge = zeros(24, 3*node_sum);
    Ge(1:3, node_id_lis(1)*3-2:node_id_lis(1)*3) = eye(3);
    Ge(4:6, node_id_lis(2)*3-2:node_id_lis(2)*3) = eye(3);
    Ge(7:9, node_id_lis(3)*3-2:node_id_lis(3)*3) = eye(3);
    Ge(10:12, node_id_lis(4)*3-2:node_id_lis(4)*3) = eye(3);
    Ge(13:15, node_id_lis(5)*3-2:node_id_lis(5)*3) = eye(3);
    Ge(16:18, node_id_lis(6)*3-2:node_id_lis(6)*3) = eye(3);
    Ge(19:21, node_id_lis(7)*3-2:node_id_lis(7)*3) = eye(3);
    Ge(22:24, node_id_lis(8)*3-2:node_id_lis(8)*3) = eye(3);
    res = Ge'*unit_stiffness;
end

