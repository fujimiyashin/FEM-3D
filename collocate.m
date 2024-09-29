lis_size = size(cubic_lis);
cubic_num = lis_size(1);
%刚度矩阵
for i=1:cubic_num
    X = [];
    Y = [];
    Z = [];
    node_id_lis = [];
    cubic = cubic_lis(i, :);
    for node=cubic
        X = [X;node.x];
        Y = [Y;node.y];
        Z = [Z;node.z];
        node_id_lis = [node_id_lis, node.id];
    end
    unit_stiffness_k = zeros(24,24);
    unit_stiffness_k = unit_stiffness_K(G, nu, X, Y, Z);
    stiffness_k = stiffness_k + total_stiffness_K(unit_stiffness_k, node_sum, node_id_lis);
end

%等效节点体积力
for i=1:cubic_num
    force = [0;0;0];
    X = [];
    Y = [];
    Z = [];
    node_id_lis = [];
    cubic = cubic_lis(i, :);
    for node=cubic
        X = [X;node.x];
        Y = [Y;node.y];
        Z = [Z;node.z];
        node_id_lis = [node_id_lis, node.id];
    end
    unit_stiffness_p = zeros(24,1);
    unit_stiffness_p = unit_stiffness_F(X, Y, Z, force);
    stiffness_p = stiffness_p + total_stiffness_P(unit_stiffness_p, node_sum, node_id_lis);
end

%等效节点面积力
for i=1:cubic_num
    const_pnt = 1;
    gauss_pnt = [-1/sqrt(3), 1/sqrt(3)];
    weight = [1,1];
    X = [];
    Y = [];
    Z = [];
    node_id_lis = [];
    cubic = cubic_lis(i, :);
    for node=cubic
        X = [X;node.x];
        Y = [Y;node.y];
        Z = [Z;node.z];
        node_id_lis = [node_id_lis, node.id];
    end

    if node_id_lis(1) == 28 || node_id_lis(1) == 29 || node_id_lis(1) == 30 || node_id_lis(1) == 61 || node_id_lis(1) == 62 || node_id_lis(1) == 63 
        force = magnitude*[0;0;-1];
    else
        force = [0;0;0];
    end
    
    unit_stiffness_p = zeros(24,1);
    unit_stiffness_p = unit_stiffness_T(gauss_pnt, const_pnt, gauss_pnt, weight, X, Y, Z, force);
    stiffness_p = stiffness_p + total_stiffness_P(unit_stiffness_p, node_sum, node_id_lis);
end