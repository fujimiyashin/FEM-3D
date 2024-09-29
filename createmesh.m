for k=1:len_z+1
    for j=1:len_y+1
        for i=1:len_x+1
            id = i+(len_x+1)*(j-1)+(len_x+1)*(len_y+1)*(k-1);
            node = Node(id, i, j, k);
            node_lis = [node_lis, node];
        end
    end
end