run init.m;
run createmesh.m;
run createcubic.m;
%run plot_cubic.m;
run collocate.m;
run apply_boundary.m;
run solve.m;

for i=1:node_sum
    x = node_lis(i).x+qe(i*3-2);
    y = node_lis(i).y+qe(i*3-1);
    z = node_lis(i).z+qe(i*3);
    node_lis(i) = Node(i, x, y, z);
end
run createcubic.m;
run plot_cubic.m;