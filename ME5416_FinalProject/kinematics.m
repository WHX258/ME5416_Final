%% Function handlers
% Rotation matrix from D-H table
A = @(a, alpha, d, q) ...
        ([cos(q) -sin(q)*cos(alpha)  sin(q)*sin(alpha) a*cos(q); ...
          sin(q)  cos(q)*cos(alpha) -cos(q)*sin(alpha) a*sin(q); ...
          0       sin(alpha)         cos(alpha)        d; ...
          0       0                  0                 1]);
% Added mass of cylindrical body
Addox = @(m, r, l) (diag([0.1*m; pi*rho*r^2*l; pi*rho*r^2*l; 0; (pi*rho*r^2*l^3)/12; (pi*rho*r^2*l^3)/12]));
Addoy = @(m, r, l) (diag([pi*rho*r^2*l; 0.1*m; pi*rho*r^2*l; (pi*rho*r^2*l^3)/12; 0; (pi*rho*r^2*l^3)/12]));
Addoz = @(m, r, l) (diag([pi*rho*r^2*l; pi*rho*r^2*l; 0.1*m; (pi*rho*r^2*l^3)/12; (pi*rho*r^2*l^3)/12; 0]));

% Rotation matrix
Rot = @(eta) ...
        [cos(eta(6))*cos(eta(5)) sin(eta(6))*cos(eta(5)) -sin(eta(5)); ...
        -sin(eta(6))*cos(eta(4)) + cos(eta(6))*sin(eta(5))*sin(eta(4)) ...
        cos(eta(6))*cos(eta(4)) + sin(eta(6))*sin(eta(5))*sin(eta(4)) ...
        sin(eta(4))*cos(eta(5));
        sin(eta(6))*sin(eta(4)) + cos(eta(6))*sin(eta(5))*cos(eta(4)) ...
        -cos(eta(6))*sin(eta(4)) + sin(eta(6))*sin(eta(5))*cos(eta(4)) ...
        cos(eta(4))*cos(eta(5))];
% Inertia tensors of cylindrical body
Iox = @(m, r, l) (m.*diag([0.5*r^2, (3*r^2 + l^2)/12, (3*r^2 + l^2)/12]));
Ioy = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, 0.5*r^2, (3*r^2 + l^2)/12]));
Ioz = @(m, r, l) (m.*diag([(3*r^2 + l^2)/12, (3*r^2 + l^2)/12, 0.5*r^2]));

%% Homogeneous transformations solution
Tr = cell(n,1);
Tr0 = eye(4,'sym'); Tr0(1:3,1:3) = Rot0 * Rot(eta); Tr0(1:3,4) = eta(1:3);
T = Tr0;
for i = 1:n
    T = T * A(DH(i,1),DH(i,2),DH(i,3),DH(i,4));
    Tr{i} = simplify(T);
end

%% Mass centers for each link
r_c_m = cell(n,1);
for i = 1:n    
    temp = Tr{i}*[[1;0;0;0] [0;1;0;0] [0;0;1;0] [c{i};1]];
	r_c_m{i} = temp(1:3,4);  
end

%% Inertia tensor/Added mass for each link relative to the inertial frame 
% computation depends on COM offset
I = cell(1,n);
A = cell(1,n);
for i = 1:n
    [~, index] = max(abs(c{i}));
    if (index == 1) 
        I{i} = Iox(m(i), r(i), l(i));
        if (rho == 0)
            A{i} = zeros(6);
        else
            A{i} = Addox(m(i), r(i), l(i));
        end
    elseif (index == 2)
        I{i} = Ioy(m(i), r(i), l(i));
        if (rho == 0)
            A{i} = zeros(6);
        else
            A{i} = Addoy(m(i), r(i), l(i));
        end
    else
        I{i} = Ioz(m(i), r(i), l(i));
        if (rho == 0)
            A{i} = zeros(6);
        else    
            A{i} = Addoz(m(i), r(i), l(i));
        end
    end
end

%% Velocity Jacobians
Jv = cell(1,n);
Jw = cell(1,n);
for i = 1:n
    z = Tr0(1:3,3);
    o = Tr0(1:3,4);
    Jvt = sym(zeros(3,n));
    Jwt = sym(zeros(3,n));
    for j = 1:i
        if (sigma(i) == 0)
            Jvt(:,j) = cross(z, r_c_m{i}-o);
            Jwt(:,j) = z;
        else
            Jvt(:,j) = z;
            Jwt(:,j) = 0;
        end
        z = Tr{j}(1:3,3);
        o = Tr{j}(1:3,4);
    end
    Jv{i} = simplify(Jvt);
    Jw{i} = simplify(Jwt);
end

Jvt_ee = sym(zeros(3,n));
for i=1:n
    Jvt_ee(:,i) = diff(Tr{n}(1:3,4), q(i));
end
Jv_ee = simplify(Jvt_ee);
Jw_ee = Jw{n}; 

%% Potential energy solution / Inertia matrix solution
PE = 0;
M_sym = 0;
for i = 1:n
    R = Tr{i}(1:3,1:3);
    M_sym = M_sym + (Jv{i}'*(m(i)*eye(3)+A{i}(1:3,1:3))*Jv{i} + ...
        Jw{i}'*R*(I{i}+A{i}(4:6,4:6))*R'*Jw{i});
    PE = PE + (m(i)*g - B(i))*r_c_m{i}(3);
end

%% The Christoffel symbols
c = zeros(n,n,n,'sym');
for k = 1:n
    for i = 1:n
        for j = 1:n
            c(i,j,k) = 0.5 * (diff(M_sym(k,j),q(i)) + ...
                diff(M_sym(k,i),q(j)) - diff(M_sym(i,j),q(k)));
        end
    end
end

%% The Coriolis matrix
C_sym = zeros(n,n,'sym');
for k = 1:n
    for j = 1:n
        temp = 0;
        for i = 1:n
            temp = temp + c(i,j,k)*dq(i);
        end
        C_sym(k,j) = temp;
    end
end

%% Friction terms
D_sym = mu.*eye(n,'sym');

%% The gravitation terms
g_sym = zeros(n,1,'sym');
for k = 1:n
    g_sym(k) = diff(PE,q(k));
end

%% Generate numeric functions
matlabFunction(M_sym,'File','param\get_M','Vars',{q});
matlabFunction(C_sym,'File','param\get_C','Vars',{[q;dq]});
matlabFunction(D_sym,'File','param\get_D','Vars',{[q;dq]});
matlabFunction(g_sym,'File','param\get_g','Vars',{q});

matlabFunction(Jv_ee,'File','param\get_Jv_ee','Vars',{q});
matlabFunction(Jw_ee,'File','param\get_Jw_ee','Vars',{q});
matlabFunction([Jv_ee;Jw_ee],'File','param\get_J_ee','Vars',{q});

matlabFunction(Tr{n}(1:3,4),'File','param\get_X_ee','Vars',{q});
