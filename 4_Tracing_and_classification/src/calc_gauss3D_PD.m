function y = calc_gauss3D_PD(x, xdata)
%Calculate 3D gaussian f(v) = exp(-v'*A*v)
%based on code by R.Xu, UCLA, 2014

    [L, M, N] = size(xdata.x); Num = L*M*N;
    v = [reshape(xdata.x - x(3),1,Num); reshape(xdata.y - x(4),1,Num); reshape(xdata.z - x(5),1,Num)];

    vector1 = [0 0 1];
    rotmat1 = MatrixQuaternionRot(vector1,x(9));
    
    vector2 = [0 1 0];
    rotmat2 = MatrixQuaternionRot(vector2,x(10));

    vector3 = [0 0 1];
    rotmat3 = MatrixQuaternionRot(vector3,x(11));

    rotMAT =  rotmat3*rotmat2*rotmat1;

    D = [x(6)  0   0;
         0  x(7)  0;
         0    0  x(8)];

    A = rotMAT' * D * rotMAT; 
    y = x(2)*reshape(exp(-dot(v,A*v,1)),L,M,N) + x(1);    
    
end



function dd = MatrixQuaternionRot(vector,theta)

theta = theta*pi/180;
vector = vector/sqrt(dot(vector,vector));
w = cos(theta/2); x = -sin(theta/2)*vector(1); y = -sin(theta/2)*vector(2); z = -sin(theta/2)*vector(3);
RotM = [1-2*y^2-2*z^2 2*x*y+2*w*z 2*x*z-2*w*y;
      2*x*y-2*w*z 1-2*x^2-2*z^2 2*y*z+2*w*x;
      2*x*z+2*w*y 2*y*z-2*w*x 1-2*x^2-2*y^2;];

dd = RotM;
end