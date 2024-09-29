fixed_node = [1, 2, 3, 34, 35, 36, 67, 68, 69];

for i=fixed_node
    stiffness_k(i*3-2:i*3, :) = 0;
    stiffness_k(:, i*3-2:i*3) = 0;
    stiffness_k(i*3-2:i*3, i*3-2:i*3) = eye(3);
    stiffness_p(i*3-2:i*3) = 0;
end