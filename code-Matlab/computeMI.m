function mi = computeMI(x,y)
joint_entropy = entropy([x;y]);

% Compute the individual entropies of x and y
entropy_x = entropy(x);
entropy_y = entropy(y);

% Compute the mutual information between x and y
mi = entropy_x + entropy_y - joint_entropy;


