function C_sym = get_C(in1)
%get_C
%    C_sym = get_C(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    18-Apr-2025 14:36:53

dq1 = in1(3,:);
dq2 = in1(4,:);
q1 = in1(1,:);
q2 = in1(2,:);
t2 = cos(q1);
t3 = sin(q1);
t4 = q1+q2;
t5 = cos(t4);
t6 = sin(t4);
t7 = t5./2.0;
t8 = t6./2.0;
t11 = dq2.*t5.*t6.*pi.*2.499325;
t9 = t2+t7;
t10 = t3+t8;
t12 = -t11;
t13 = t5.*t10.*pi.*1.485e-2;
t14 = t6.*t9.*pi.*5.0135;
t15 = -t13;
t17 = -dq2.*(t13-t14);
t18 = dq2.*(t13-t14);
t16 = t14+t15;
C_sym = reshape([t18-dq1.*(t2.*t3.*pi.*2.499325+t9.*t10.*pi.*9.9973),t12-dq1.*(t5.*t10.*pi.*5.0135-t6.*t9.*pi.*1.485e-2),t18+dq1.*(t13-t14),t12-dq1.*t5.*t6.*pi.*2.499325],[2,2]);
