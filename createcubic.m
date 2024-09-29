cubic_lis = [];
for k=1:len_z
    for j=1:len_y
        for i=1:len_x
            id = i+(len_x+1)*(j-1)+(len_x+1)*(len_y+1)*(k-1);
            cubic = [node_lis(id), node_lis(id+1), node_lis(id+len_x+2), node_lis(id+len_x+1), node_lis(id+(len_x+1)*(len_y+1)), node_lis((len_x+1)*(len_y+1)+id+1), node_lis((len_x+1)*(len_y+1)+id+len_x+2), node_lis((len_x+1)*(len_y+1)+id+len_x+1)];
            cubic_lis = [cubic_lis;cubic];
        end
    end
end