function out1 = get_J_ee(in1)
%get_J_ee
%    OUT1 = get_J_ee(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    18-Apr-2025 14:36:53

q1 = in1(1,:);
q2 = in1(2,:);
t2 = q1+q2;
t3 = cos(t2);
t4 = sin(t2);
t5 = -t4;
out1 = reshape([t5-sin(q1),0.0,t3+cos(q1),0.0,-1.0,0.0,t5,0.0,t3,0.0,-1.0,0.0],[6,2]);
