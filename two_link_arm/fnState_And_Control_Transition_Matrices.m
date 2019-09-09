function [A,B] = fnState_And_Control_Transition_Matrices(x,u,du,dt)

global d1;
global d2;
global d3;
global b1;
global b2;
global b1_2;
global b2_1;

x1 = x(1,1);
x2 = x(2,1);
x3 = x(3,1);
x4 = x(4,1);

u1 = u(1,1);
u2 = u(2,1);


det = d1 * d3  - d3^2 - (d2 * cos(x2))^2;



A = zeros(4,4);


B = zeros(4,2);

A(1,3) = 1;
A(2,4)  = 1;


         
A(3,2) = -2/(d1*d3-d3^2-cos(x2)^2)^2*(-d2*d3*(x3+x4)^2*sin(x2)-d2^2*x3^2*sin(x2)*cos(x2)-d2*(b2_1*x3+b2*x4)*cos(x2)+(d3*b1-d3*b2_1)*x3+(d3*b1_2-d3*b2)*x4)*cos(x2)*sin(x2)+1/(d1*d3-d3^2-cos(x2)^2)*(-d2*d3*(x3+x4)^2*cos(x2)-d2^2*x3^2*cos(x2)^2+d2^2*x3^2*sin(x2)^2+d2*(b2_1*x3+b2*x4)*sin(x2));
A(3,3) = 1/(d1*d3-d3^2-cos(x2)^2)*(-2*d2*d3*(x3+x4)*sin(x2)-2*d2^2*x3*sin(x2)*cos(x2)-d2*b2_1*cos(x2)+d3*b1-d3*b2_1);
A(3,4) = 1/(d1*d3-d3^2-cos(x2)^2)*(-2*d2*d3*(x3+x4)*sin(x2)-d2*b2*cos(x2)+d3*b1_2-d3*b2); 
A(4,2) = -2/(d1*d3-d3^2-cos(x2)^2)^2*(d2*d3*x4*(2*x3+x4)*sin(x2)+d1*d2*x3^2*sin(x2)+d2^2*(x3+x4)^2*sin(x2)*cos(x2)+d2*((2*b2_1-b1)*x3+(2*b2-b1_2)*x4)*cos(x2)+(d1*b2_1-d3*b1)*x3+(d1*b2-d3*b1_2)*x4)*cos(x2)*sin(x2)+1/(d1*d3-d3^2-cos(x2)^2)*(d2*d3*x4*(2*x3+x4)*cos(x2)+d1*d2*x3^2*cos(x2)+d2^2*(x3+x4)^2*cos(x2)^2-d2^2*(x3+x4)^2*sin(x2)^2-d2*((2*b2_1-b1)*x3+(2*b2-b1_2)*x4)*sin(x2));
A(4,3) = 1/(d1*d3-d3^2-cos(x2)^2)*(2*d2*d3*x4*sin(x2)+2*d1*d2*x3*sin(x2)+2*d2^2*(x3+x4)*sin(x2)*cos(x2)+d2*(2*b2_1-b1)*cos(x2)+d1*b2_1-d3*b1);
A(4,4) = 1/(d1*d3-d3^2-cos(x2)^2)*(d2*d3*(2*x3+x4)*sin(x2)+d2*d3*x4*sin(x2)+2*d2^2*(x3+x4)*sin(x2)*cos(x2)+d2*(2*b2-b1_2)*cos(x2)+d1*b2-d3*b1_2);
 
 

B(3,1) = d3;
B(3,2) = -(d3 + d2 * cos(x2));
B(4,1) = B(3,2);
B(4,2) = d1 + 2 * d2 * cos(x2);

B = 1/det * B;
c_1 = 1/det * [ 0 0 0 d3 * u1]';
c_2 = 1/det * [0, 0, -d3 + d2 * cos(x2)*u1,d1 + 2 * d2 * cos(x2) * u2]';

C_1 = zeros(4,2);
C_1(3,1) = d3;
C_1(4,2) = d2 * cos(x2);
C_1 = C_1 * 1/det;

C_2 = zeros(4,2);
C_2(3,1) = d2 * cos(x2);
C_2(4,2) = 2* d2 * cos(x2);
C_2 = C_1 * 1/det;

C = [C_1 C_2] *sqrt(dt);
c = [c_1 c_2]*sqrt(dt);