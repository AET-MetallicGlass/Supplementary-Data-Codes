function dd = MatrixQuaternionRot(vector,theta)

theta = theta*pi/180;
% vector = vector/sqrt(dot(vector,vector));
vector = vector ./ norm(vector);
w = cos(theta/2); x = -sin(theta/2)*vector(1); y = -sin(theta/2)*vector(2); z = -sin(theta/2)*vector(3);
RotM = [1-2*y^2-2*z^2 2*x*y+2*w*z 2*x*z-2*w*y;
      2*x*y-2*w*z 1-2*x^2-2*z^2 2*y*z+2*w*x;
      2*x*z+2*w*y 2*y*z-2*w*x 1-2*x^2-2*y^2;];

dd = RotM;
end