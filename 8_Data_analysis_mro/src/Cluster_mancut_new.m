classdef Cluster_mancut_new

  properties
    origin
    N_atoms
    pos
    idxs
    type
    param
    grid_pts
    r_grid_cutoff = 1/3.8;
    N_atoms_relax_off = 7 %3300;
  end
  
  methods
    function obj = Cluster_mancut_new(origin_pos,idx,type_in,cut)
      obj.origin   = origin_pos;
      obj.N_atoms  = 1;
      obj.pos      = obj.origin;
      obj.idxs     = idx;
      obj.type     = type_in;
      obj.param    = [1,1,0,0,0];
      obj.grid_pts = zeros(1,size(nn_basis(type_in),1));
      obj.r_grid_cutoff = cut;
    end

    % Returns obj with new_pos attempted. If not successful, then code will
    % return 0, and behavior of obj is undefined.
    % If addition of new_pos is successful, then obj is an updated cluster.
    function [obj,code] = add_pos(obj,new_pos,new_idx)
      if(obj.N_atoms == 1)
        obj.N_atoms  = 2;
        obj.pos  = [obj.pos; new_pos];
        obj.idxs = [obj.idxs; new_idx];
        obj.grid_pts(2,2) = 1;
        obj  = relax_vars(obj);
        code = 1;
        return;
      end
      if(obj.type <= 4)
        [obj, code] = add_pos_crystal(obj, new_pos, new_idx);
        return;
      elseif(obj.type == 5)
        [obj, code] = add_pos_nn(obj, new_pos, new_idx);
        return;
      else
        code = 0;
        return;
      end
      return;
    end

    function plot_cluster(obj)
      r_grid = calc_r_grid(obj.param,obj.grid_pts,obj.type);
      r = obj.pos-obj.origin;
      scatter3(r(:,1), r(:,2), r(:,3),120,'magenta','filled');
      hold on;
      scatter3(r_grid(:,1), r_grid(:,2), r_grid(:,3),...
        150, 'green', 'filled' ...
      );
      % for i=2:obj.N_atoms
      %   line([r(i,1),r_grid(i,1)],[r(i,2),r_grid(i,2)],[r(i,3),r_grid(i,3)],'linewidth',6);
      % end
      adj = neigh_matrix(obj);
      for i=1:obj.N_atoms
        for j=i+1:obj.N_atoms
          if(adj(i,j)==1)
            line([r_grid(i,1),r_grid(j,1)],[r_grid(i,2),r_grid(j,2)],[r_grid(i,3),r_grid(j,3)],'linewidth',4);
            % line([r(i,1),r(j,1)],[r(i,2),r(j,2)],[r(i,3),r(j,3)],'linewidth',4);
          end
        end
      end
      axis equal;
      xlim([-4.8,4.8]);
      ylim([-4.8,4.8]);
      zlim([-4.8,4.8]);
      hold off;
    end

    function adj = neigh_matrix(obj)
      r_grid_pts = calc_r_grid(obj.param,obj.grid_pts,obj.type)/obj.param(1);
      for i=1:obj.N_atoms
        for j=i+1:obj.N_atoms
          r_ij = r_grid_pts(i,:)-r_grid_pts(j,:);
          d_ij = sum(r_ij.*r_ij,2);
          if(obj.type<=4 && abs(d_ij-1) < 0.001)
            adj(i,j) = 1;
            adj(j,i) = 1;
          elseif(obj.type==5 && abs(d_ij-1.05) < 0.1)
            adj(i,j) = 1;
            adj(j,i) = 1;
          end
        end
      end
    end

    % Returns the L2 norm of the positions with vars and position
    function err = calc_err_obj(obj)
      r_grid = calc_r_grid(obj.param,obj.grid_pts,obj.type);
      r_diff = obj.pos - obj.origin - r_grid;
      err    = sum(sum(r_diff.*r_diff,2));
    end

  end
  
end

% Relaxes the 5 variables that best matches positions at pos
function cluster = relax_vars(cluster)
  N_ITER_MAX     = 1000;
  DER_CUTOFF     = 0.01;
  GAMMA          = 0.01;
  STEP_A         = 0.01;
  STEP_THETA     = 0.01;
  NOISE          = 0.0;

  pos      = cluster.pos;
  for i=1:cluster.N_atoms
    pos(i,:) = pos(i,:)-cluster.origin;
  end
  type     = cluster.type;
  vars     = cluster.param;
  grid_pts = cluster.grid_pts;
  for n_iter=1:N_ITER_MAX
    vars_a_plus  = vars+[STEP_A,0,0,0,0];
    vars_a_minus = vars+[-1*STEP_A,0,0,0,0];
    d_err_d_a = (calc_err(vars_a_plus,pos,grid_pts,type)-calc_err(vars_a_minus,pos,grid_pts,type))/(2*STEP_A);
    del_a = -1*GAMMA*d_err_d_a+(rand*2-1)*NOISE;

    qx_plus  = axis_angle_to_quat([1,0,0],STEP_THETA);
    qx_minus = axis_angle_to_quat([1,0,0],-1*STEP_THETA);
    vars_qx_plus  = [vars(1), quat_product(qx_plus, vars(2:5))];
    vars_qx_minus = [vars(1), quat_product(qx_minus,vars(2:5))];
    d_err_d_qx = (calc_err(vars_qx_plus,pos,grid_pts,type)-calc_err(vars_qx_minus,pos,grid_pts,type))/(2*STEP_THETA);
    del_thetax = -1*GAMMA*d_err_d_qx+(rand*2-1)*NOISE;

    qy_plus  = axis_angle_to_quat([0,1,0],STEP_THETA);
    qy_minus = axis_angle_to_quat([0,1,0],-1*STEP_THETA);
    vars_qy_plus  = [vars(1), quat_product(qy_plus, vars(2:5))];
    vars_qy_minus = [vars(1), quat_product(qy_minus,vars(2:5))];
    d_err_d_qy = (calc_err(vars_qy_plus,pos,grid_pts,type)-calc_err(vars_qy_minus,pos,grid_pts,type))/(2*STEP_THETA);
    del_thetay = -1*GAMMA*d_err_d_qy+(rand*2-1)*NOISE;

    qz_plus  = axis_angle_to_quat([0,0,1],STEP_THETA);
    qz_minus = axis_angle_to_quat([0,0,1],-1*STEP_THETA);
    vars_qz_plus  = [vars(1), quat_product(qz_plus, vars(2:5))];
    vars_qz_minus = [vars(1), quat_product(qz_minus,vars(2:5))];
    d_err_d_qz = (calc_err(vars_qz_plus,pos,grid_pts,type)-calc_err(vars_qz_minus,pos,grid_pts,type))/(2*STEP_THETA);
    del_thetaz = -1*GAMMA*d_err_d_qz+(rand*2-1)*NOISE;

    q_new = quat_product(axis_angle_to_quat([1,0,0],del_thetax),vars(2:5));
    q_new = quat_product(axis_angle_to_quat([0,1,0],del_thetay),q_new);
    q_new = quat_product(axis_angle_to_quat([0,0,1],del_thetaz),q_new);
    vars = [vars(1)+del_a, q_new];
    der_tot = sqrt( ...
      d_err_d_a*d_err_d_a + ...
      d_err_d_qx*d_err_d_qx + d_err_d_qy*d_err_d_qy + d_err_d_qz*d_err_d_qz ...
    );

    if(der_tot < DER_CUTOFF)
      break;
    end
  end
  cluster.param = vars;
end

function [obj, code] = add_pos_nn(obj, new_pos, new_idx)
  basis_N = size(nn_basis(obj.type),1);
  r       = new_pos-obj.origin;

  % Determine new grid point
  new_grid_pt = zeros(1,size(nn_basis(obj.type),1));
%   found_nearest_grid_pt = 0;
  grid_trials = repmat(new_grid_pt,basis_N,1);
  grid_trials = grid_trials + eye(basis_N);
  r_trials = calc_r_grid(obj.param,grid_trials,obj.type);
  for i=1:basis_N
    r_trials(i,:) = r_trials(i,:)-r;
  end
  [~,idx] = min(sum(r_trials.*r_trials,2));
  new_grid_pt = grid_trials(idx,:);

  % Determine continuity for new grid point. If not continous, then return
  % code as 0
  r_new_grid_pt  = calc_r_grid(obj.param,new_grid_pt,obj.type);
  r_grid_pts     = calc_r_grid(obj.param,obj.grid_pts,obj.type);
  r_diff  = r_grid_pts - repmat(r_new_grid_pt,size(r_grid_pts,1),1);
  r_diff_mag = sqrt(sum(r_diff.*r_diff,2))/obj.param(1);
  min_r_diff_mag = min(r_diff_mag);
  if(abs(min_r_diff_mag-1) < 0.001)
    obj.N_atoms  = obj.N_atoms + 1;
    obj.pos      = [obj.pos; new_pos];
    obj.idxs     = [obj.idxs; new_idx];
    obj.grid_pts = [obj.grid_pts; new_grid_pt];
  else
    code = 0;
    return;
  end

  obj = relax_vars(obj);
  r_grid = calc_r_grid(obj.param,obj.grid_pts,obj.type);
%   for i=1:size(r_diff,1)
%     r_diff(i,:) = obj.pos(i,:) - obj.origin - r_grid(i,:);
%   end
  r_diff = obj.pos - repmat(obj.origin, [size(obj.pos,1),1]) - r_grid;
  dists  = sqrt(sum(r_diff.*r_diff,2));
  if(max(dists)>obj.r_grid_cutoff)
    code = 0;
  else
    code = 1;
  end
  return;
end

function [obj, code] = add_pos_crystal(obj, new_pos, new_idx)
  basis_N = size(nn_basis(obj.type),1);
  r       = new_pos-obj.origin;
  
  % Determine new grid point
  new_grid_pt = zeros(1,size(nn_basis(obj.type),1));
  found_nearest_grid_pt = 0;
  while(~found_nearest_grid_pt)
    grid_trials = repmat(new_grid_pt,basis_N,1);
    grid_trials = grid_trials + eye(basis_N);
    r_trials = calc_r_grid(obj.param,grid_trials,obj.type);
    for i=1:basis_N
      r_trials(i,:) = r_trials(i,:)-r;
    end
    [~,idx] = min(sum(r_trials.*r_trials,2));
    if(idx==1)
      found_nearest_grid_pt = 1;
    end
    new_grid_pt = grid_trials(idx,:);
  end

  % Determine continuity for new grid point. If not continous, then return
  % code as 0
  r_new_grid_pt  = calc_r_grid(obj.param,new_grid_pt,obj.type);
  r_grid_pts     = calc_r_grid(obj.param,obj.grid_pts,obj.type);
  r_diff  = r_grid_pts - repmat(r_new_grid_pt,size(r_grid_pts,1),1);
  r_diff_mag = sqrt(sum(r_diff.*r_diff,2))/obj.param(1);
  min_r_diff_mag = min(r_diff_mag);
  if(abs(min_r_diff_mag-1) < 0.001)
    obj.N_atoms  = obj.N_atoms + 1;
    obj.pos      = [obj.pos; new_pos];
    obj.idxs     = [obj.idxs; new_idx];
    obj.grid_pts = [obj.grid_pts; new_grid_pt];
  else
    code = 0;
    return;
  end

  % If less atoms than N_atoms_relax_off, then check if every atom is below
  % cutoff after relaxing variables
  if(obj.N_atoms <= obj.N_atoms_relax_off)
    obj = relax_vars(obj);
  end
    r_grid = calc_r_grid(obj.param,obj.grid_pts,obj.type);
%     for i=1:size(r_diff,1)
%       r_diff(i,:) = obj.pos(i,:) - obj.origin - r_grid(i,:);
%     end
    r_diff = obj.pos - repmat(obj.origin, [size(obj.pos,1),1]) - r_grid;
    dists  = sqrt(sum(r_diff.*r_diff,2));
    if(max(dists)>obj.r_grid_cutoff)
      code = 0;
    else
      code = 1;
    end
  return;
end

% Returns the L2 norm of the positions with vars and position
function err = calc_err(vars,pos,grid_pts,type)
  r_grid = calc_r_grid(vars,grid_pts,type);
  r_diff = pos - r_grid;
  err    = sum(sum(r_diff.*r_diff,2));
end

function r_grid = calc_r_grid(vars,grid_pts,type)
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

function basis = nn_basis(type)
  if(type==1)
    N_grid_pts = 7;
    basis = [ ...
      0,0,0  ;...
      1,0,0  ;...
      -1,0,0 ;...
      0,1,0  ;...
      0,-1,0 ;...
      0,0,1  ;...
      0,0,-1 ;...
    ];
  elseif(type==2)
    d = 1/sqrt(3);
    basis = [ ...
      0,0,0  ;...
      d,d,d  ;...
      -d,d,d ;...
      d,-d,d ;...
      d,d,-d ;...
      -d,-d,d ;...
      -d,d,-d ;...
      d,-d,-d ;...
      -d,-d,-d ;...
    ];
  elseif(type==3)
    d = 1/sqrt(2);
    basis = [ ...
      0,0,0    ;...
      d,d,0    ;...
      -d,d,0   ;...
      d,-d,0   ;...
      -d,-d,0  ;...
      d,0,d    ;...
      -d,0,d   ;...
      d,0,-d   ;...
      -d,0,-d  ;...
      0,d,d    ;...
      0,-d,d   ;...
      0,d,-d   ;...
      0,-d,-d  ;...
    ];
  elseif(type==4)
    basis = [ ...
      0,0,0                           ;...
      1/2,1/(2*sqrt(3)),sqrt(2/3)     ;...
      -1/2,1/(2*sqrt(3)),sqrt(2/3)    ;...
      0,-1/sqrt(3),sqrt(2/3)          ;...
      1/2,1/(2*sqrt(3)),-sqrt(2/3)    ;...
      -1/2,1/(2*sqrt(3)),-sqrt(2/3)   ;...
      0,-1/sqrt(3),-sqrt(2/3)         ;...
      1,0,0                           ;...
      1/2,sqrt(3)/2,0                 ;...
      -1/2,sqrt(3)/2,0                ;...
      -1,0,0                          ;...
      -1/2,-sqrt(3)/2,0               ;...
      1/2,-sqrt(3)/2,0                ;...
    ];
  elseif(type==5)
    phi = (1+sqrt(5))/2;
    den = sqrt(1+phi*phi);
    basis = [ ...
      0,0,0       ;...
      0,1,phi     ;...
      0,-1,phi    ;...
      0,1,-phi    ;...
      0,-1,-phi   ;...
      1,phi,0     ;...
      -1,phi,0    ;...
      1,-phi,0    ;...
      -1,-phi,0   ;...
      phi,0,1     ;...
      -phi,0,1    ;...
      phi,0,-1    ;...
      -phi,0,-1   ;...
    ]/den;
  else
    basis = [];
  end
end

function m = quat_to_rotm(q_raw)
  q = q_raw/norm(q_raw);
  m = [ ...
    1-2*(q(3)*q(3)+q(4)*q(4)), 2*(q(2)*q(3)-q(4)*q(1)), 2*(q(2)*q(4)+q(3)*q(1))  ; ...
    2*(q(2)*q(3)+q(4)*q(1)), 1-2*(q(2)*q(2)+q(4)*q(4)), 2*(q(3)*q(4)-q(2)*q(1))  ; ...
    2*(q(2)*q(4)-q(3)*q(1)), 2*(q(3)*q(4)+q(2)*q(1)), 1-2*(q(2)*q(2)+q(3)*q(3)) ...
  ];
end

function q = axis_angle_to_quat(axis_raw,angle)
  axis = axis_raw/norm(axis_raw);
  q = [cos(angle/2),axis*sin(angle/2)];
end

% Returns r = p*q where p and q are quaternions
function r = quat_product(p_raw,q_raw)
  p = p_raw/norm(p_raw);
  q = q_raw/norm(q_raw);
  r = zeros(4,1);
  p0    = p(1);    q0    = q(1);
  p_vec = p(2:4);  q_vec = q(2:4);
  r(1) = p0*q0-dot(p_vec,q_vec);
  r(2:4) = p0*q_vec+q0*p_vec+cross(p_vec,q_vec);
  r = r';
end

