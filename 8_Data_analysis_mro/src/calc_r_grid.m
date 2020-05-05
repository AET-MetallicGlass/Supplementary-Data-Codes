function [r_grid,v] = calc_r_grid(vars,grid_pts,type)
  if(type ~= 4)
    a = vars(1);
    q = vars(2:5);
    rotm = quat_to_rotm(q);
    v = a*(nn_basis(type)')*(grid_pts');
    r_grid = (rotm*v)';
  else
    basis = nn_basis(type);
    v = zeros(size(grid_pts,1),3);
    for i=1:size(grid_pts,1)
      curr_grid_pt = grid_pts(i,:); 
      level = 0;
      while(sum(curr_grid_pt)>0)
        j = 1;
        while(curr_grid_pt(j) == 0)
          j = j+1;
        end
        if(j==2 || j==3 || j==4)
          b = basis(j,:);
          b(2) = (-1)^(level)*b(2);
          v(i,:) = v(i,:) + b;
          level = level + 1;
        elseif(j==5 || j==6 || j==7)
          b = basis(j,:);
          b(2) = (-1)^(level)*b(2);
          v(i,:) = v(i,:) + b;
          level = level - 1;
        else
          v(i,:) = v(i,:) + basis(j,:);
        end
        curr_grid_pt(j) = curr_grid_pt(j)-1;
      end
    end
    a = vars(1);
    q = vars(2:5);
    rotm = quat_to_rotm(q);
    v = a*v';
    r_grid = (rotm*v)';
  end
end