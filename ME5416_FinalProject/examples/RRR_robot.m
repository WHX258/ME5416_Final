n = 3; % link number
sigma = [0; 0; 0]; % 0 if revolute, 1 if prismatic
mu = [0.5; 0.5; 0.5]; % dry friction
l = [0; 1; 1]; % link length [m]
r = [0.1; 0.1; 0.1]; % radius [m]
m = [0; 1; 1]; % mass of each link [kg]
B = [0; 0; 0]; % buoyancy [N] 
Rot0 = eye(3); % rotation matrix

eta = [ 0; 0; 0; 0; 0; 0 ]; % base pose

q = sym('q', [n 1], 'real'); % q1: shoulder yaw (Z), q2: shoulder pitch (Y), q3: elbow pitch (Y)

%%%%%% a     alpha     d      theta %%%%%%
DH = [ 0      pi/2     0      q(1);  % yaw: z-rotation
       1      0        0      q(2);  % pitch: y-rotation
       1      0        0      q(3)]; % elbow: y-rotation, link2 to link3

% 连杆质心位置，相对于本地坐标系
c{1} = [-l(1)/2; 0; 0];
c{2} = [-l(2)/2; 0; 0];
c{3} = [-l(3)/2; 0; 0];

dq = sym('dq', [n 1], 'real'); % joint velocities

%% Constants
g = 9.81;
rho = 1000;

%% Initial Conditions
q0 = deg2rad([0; 45; 90]); % initial position
dq0 = zeros(n,1);