Im = im2double (imread('Input image'));
I = rgb2gray(Im);

% Canny edge detection with emprirical threshold
edgs = edge(I, 'Canny', 0.06 );
figure(1), imshow([edgs])

% create the hough transform from our Binary image
[H, T, R] = hough(edgs);
figure(2) 
imshow(imadjust(H),[], ...
           'XData',T, ... 
           'YData',R,...
           'InitialMagnification','fit');
xlabel('\theta(degrees)') 
ylabel('\rho')
axis on 
axis normal 
hold on
colormap(gca, hot)

% Find peaks in the Hough Transform 
P = houghpeaks(H, 70, 'threshold', ceil(0.25*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
plot(x, y, 's', 'color', 'black');

% Find lines and plot them
lines = houghlines(edgs, T, R, P, 'FillGap', 5, 'MinLength', 7);
figure(3), imshow(Im), hold on
max_len = 0;
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:,1), xy(:,2), 'LineWidth',2, 'Color','green');
    
    plot(xy(1,1), xy(1,2),'x', 'LineWidth', 1,'Color', 'yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',1,'Color','red');
end 

% expect 20 points 
% 1-4 points: vertical lines
% 5-8 points: horizontal lines
% 9-14 points: left side bandeon
% 15-20 points: right side bandeon 
[x, y] = getpts;

plot(x,y,'or','MarkerSize',12);

FNT_SZ = 20;
% vertical Floor lines
vert1=[x(1) y(1) 1]';
text(vert1(1), vert1(2), 'v1', 'FontSize', FNT_SZ, 'Color', 'k')
vert2=[x(2) y(2) 1]';
text(vert2(1), vert2(2), 'v2', 'FontSize', FNT_SZ, 'Color', 'k')
vert3=[x(3) y(3) 1]';
text(vert3(1), vert3(2), 'v3', 'FontSize', FNT_SZ, 'Color', 'k')
vert4=[x(4) y(4) 1]';
text(vert4(1), vert4(2), 'v4', 'FontSize', FNT_SZ, 'Color', 'k')

myline=[vert1';vert2'];
line(myline(:,1),myline(:,2),'LineWidth',5);
myline=[vert3';vert4'];
line(myline(:,1),myline(:,2),'LineWidth',5);

%horizontal floor lines
hor1=[x(5) y(5) 1]';
text(hor1(1), hor1(2), 'h1', 'FontSize', FNT_SZ, 'Color', 'k')
hor2=[x(6) y(6) 1]';
text(hor2(1), hor2(2), 'h2', 'FontSize', FNT_SZ, 'Color', 'k')
hor3=[x(7) y(7) 1]';
text(hor3(1), hor3(2), 'h3', 'FontSize', FNT_SZ, 'Color', 'k')
hor4=[x(8) y(8) 1]';
text(hor4(1), hor4(2), 'h4', 'FontSize', FNT_SZ, 'Color', 'k')

myline=[hor1';hor2'];
line(myline(:,1),myline(:,2),'LineWidth',5, 'Color', 'g');
myline=[hor3';hor4'];
line(myline(:,1),myline(:,2),'LineWidth',5, 'Color', 'g');

% convention: clockwise from front-left (a,b,c,d)
%left side of bandoneon
a=[x(9) y(9) 1]';
text(a(1), a(2), 'a', 'FontSize', FNT_SZ, 'Color', 'w')
b=[x(10) y(10) 1]';
text(b(1), b(2), 'b', 'FontSize', FNT_SZ, 'Color', 'w')
c=[x(11) y(11) 1]';
text(c(1), c(2), 'c', 'FontSize', FNT_SZ, 'Color', 'w')
d=[x(12) y(12) 1]';
text(d(1), d(2), 'd', 'FontSize', FNT_SZ, 'Color', 'w')
e=[x(13) y(13) 1]';
text(e(1), e(2), 'e', 'FontSize', FNT_SZ, 'Color', 'w')
h=[x(14) y(14) 1]';
text(h(1), h(2), 'h', 'FontSize', FNT_SZ, 'Color', 'w')

myline=[a';b';c';d';a'];
line(myline(:,1),myline(:,2),'LineWidth',5, 'Color', 'r');
myline=[a';e';h';d'];
line(myline(:,1),myline(:,2),'LineWidth',5, 'Color', 'r');

%right side of bandoneon
a2=[x(15) y(15) 1]';
text(a2(1), a2(2), 'a2', 'FontSize', FNT_SZ, 'Color', 'w')
b2=[x(16) y(16) 1]';
text(b2(1), b2(2), 'b2', 'FontSize', FNT_SZ, 'Color', 'w')
c2=[x(17) y(17) 1]';
text(c2(1), c2(2), 'c2', 'FontSize', FNT_SZ, 'Color', 'w')
d2=[x(18) y(18) 1]';
text(d2(1), d2(2), 'd2', 'FontSize', FNT_SZ, 'Color', 'w')
e2=[x(19) y(19) 1]';
text(e2(1), e2(2), 'e2', 'FontSize', FNT_SZ, 'Color', 'w')
h2=[x(20) y(20) 1]';
text(h2(1), h2(2), 'h2', 'FontSize', FNT_SZ, 'Color', 'w')

myline=[a2';b2';c2';d2';a2'];
line(myline(:,1),myline(:,2),'LineWidth',5, 'Color', 'r');
myline=[a2';e2';h2';d2'];
line(myline(:,1),myline(:,2),'LineWidth',5, 'Color', 'r');

% let's remove the projective distortion
% vanishing points
vp(3, :) = cross(cross(hor1, hor2), cross(hor3, hor4));
vp(4, :) = cross(cross(vert1, vert2), cross(vert3, vert4));
vp(1, :) = cross(cross(b2, a2), cross(c2, d2));
vp(2, :) = cross(cross(b, a), cross(c, d));
for i = 1:length(vp(:, 1))
    vp(i, :) = vp(i, :)/vp(i,3);
end

% according to Bob Collin's algorithm 
M = zeros(3,3);
for i = 1: length(vp(:,1))
    M(1,1) = M(1,1)+vp(i,1)*vp(i,1);
    M(1,2) = M(1,2)+vp(i,1)*vp(i,2);
    M(1,3) = M(1,3)+vp(i,1)*vp(i,3);
    M(2,1) = M(2,1)+vp(i,1)*vp(i,2);
    M(2,2) = M(2,2)+vp(i,2)*vp(i,2);
    M(2,3) = M(2,3)+vp(i,2)*vp(i,3);
    M(3,1) = M(3,1)+vp(i,1)*vp(i,3);
    M(3,2) = M(3,2)+vp(i,2)*vp(i,3);
    M(3,3) = M(3,3)+vp(i,3)*vp(i,3); 
end
[vl,~] = eigs(M, 1, 'SM');
vl = vl/vl(3);

H1 = [1, 0, 0; 0, 1, 0; vl'];
%Hinv = inv(H1); 

T = maketform('projective', H1');
CI = imtransform(Im, T, 'XYScale', 1);

figure; imshow(CI);

% remove the projective distortion for the used coordinates 
[a_proj] = tformfwd(T, a(1), a(2));
a_proj(:,3) = 1;
[b_proj] = tformfwd(T, b(1), b(2));
b_proj(:,3) = 1;
[c_proj] = tformfwd(T, c(1), c(2));
c_proj(:,3) = 1;
[d_proj] = tformfwd(T, d(1), d(2));
d_proj(:,3) = 1;

[a2_proj] = tformfwd(T, a2(1), a2(2));
a2_proj(:,3) = 1;
[b2_proj] = tformfwd(T, b2(1), b2(2));
b2_proj(:,3) = 1;
[c2_proj] = tformfwd(T, c2(1), c2(2));
c2_proj(:,3) = 1;
[d2_proj] = tformfwd(T, d2(1), d2(2));
d2_proj(:,3) = 1;

% start removing affine distortion 
l1 = cross(a_proj, b_proj);
l2 = cross(a_proj, d_proj);
l3 = cross(a2_proj, b2_proj);
l4 = cross(a2_proj, d2_proj);

% defining the variables to find S
%B = [];
B = [l1(1)*l3(1), (l1(2)*l3(1)+l1(1)*l3(2)); 
    l2(1)*l4(1), (l2(2)*l4(1)+l2(1)*l4(2));];
C = [-(l1(2)*l3(2)); -(l2(2)*l4(2))];

%solving for the linear constraint
Un = B\C;
%V  = B*Un - C;

S = [Un(1) Un(2); Un(2), 1];
[U, D, V] = svd(S);
Ds = sqrt(D);

A = V*Ds*V';
H2 = [A(1,1) A(1,2) 0;A(2,1) A(2,2) 0; 0 0 1];

if H2(1,1) < 0
    H2(1,1) = -H2(1,1);
elseif H2(2,2) < 0
    H2(2,2) = -H2(2,2);
end

Hinv = H2/H1;
H = inv(Hinv);

T2 = maketform('projective', H');
CI = imtransform(Im, T2, 'XYScale', 1);
figure; imshow(CI);

a_rect = tformfwd(T, a(1), a(2));
b_rect = tformfwd(T, b(1), b(2));
c_rect = tformfwd(T, c(1), c(2));
d_rect = tformfwd(T, d(1), d(2));
e_rect = tformfwd(T, e(1), e(2));
h_rect = tformfwd(T, e(1), e(2));
a2_rect = tformfwd(T, a2(1), a2(2));
b2_rect = tformfwd(T, b2(1), b2(2));
c2_rect = tformfwd(T, c2(1), c2(2));
d2_rect = tformfwd(T, d2(1), d2(2));
e2_rect = tformfwd(T, e2(1), e2(2));
h2_rect = tformfwd(T, h2(1), h2(2));
a_rect(:,3) = 1;
b_rect(:,3) = 1;
c_rect(:,3) = 1;
d_rect(:,3) = 1;
e_rect(:,3) = 1;
h_rect(:,3) = 1;
a2_rect(:,3) = 1;
b2_rect(:,3) = 1;
c2_rect(:,3) = 1;
d2_rect(:,3) = 1;
e2_rect(:,3) = 1;
h2_rect(:,3) = 1;
% let's start calculating omega 
%zero-skew 
v5 = [0 1 0 0 0 0];

% first using the two circular points (rectfied plane constraint)
i = 1; j = 2;
v12 = [H(1, i)*H(1,j), H(1,i)*H(2,j)+H(2,i)*H(1,j), H(2,i)*H(2,j), H(3,i)*H(1,j)+H(1,i)*H(3,j), H(3,i)*H(2,j)+H(2,i)*H(3,j), H(3,i)*H(3,j)];
i = 1; j = 1;
v11 = [H(1, i)*H(1,j), H(1,i)*H(2,j)+H(2,i)*H(1,j), H(2,i)*H(2,j), H(3,i)*H(1,j)+H(1,i)*H(3,j), H(3,i)*H(2,j)+H(2,i)*H(3,j), H(3,i)*H(3,j)];
i = 2; j = 2;
v22 = [H(1, i)*H(1,j), H(1,i)*H(2,j)+H(2,i)*H(1,j), H(2,i)*H(2,j), H(3,i)*H(1,j)+H(1,i)*H(3,j), H(3,i)*H(2,j)+H(2,i)*H(3,j), H(3,i)*H(3,j)];

% line at infinity, vanishing point relation 
v = cross(cross(a, e), cross(d, h)); 
v = v/v(3);
%v3 = [P(1)*v(1), P(1)*v(2)+P(2)*v(1), P(2)*v(2), P(1)*v(3)+P(3)*v(1), P(2)*v(3)+P(3)*v(2), P(3)*v(3)];
%v4 = [Q(1)*v(1), Q(1)*v(2)+Q(2)*v(1), Q(2)*v(2), Q(1)*v(3)+Q(3)*v(1), Q(2)*v(3)+Q(3)*v(2), Q(3)*v(3)];
v3 = [H(1, 1)*v(1), H(1, 1)*v(2)+H(2, 1)*v(1), H(2, 1)*v(2), H(1, 1)*v(3)+H(3, 1)*v(1), H(2, 1)*v(3)+H(3, 1)*v(2), H(3, 1)*v(3)];
v4 = [H(1, 2)*v(1), H(1, 2)*v(2)+H(2, 2)*v(1), H(2, 2)*v(2), H(1, 2)*v(3)+H(3, 2)*v(1), H(2, 2)*v(3)+H(3, 2)*v(2), H(3, 2)*v(3)];

V = [v5; v12; (v11-v22); v3; v4;];
[U, D, T] = svd(V); 
omega  = T(:,6);

w =[omega(1) omega(2) omega(4); omega(2) omega(3) omega(5); omega(4) omega(5) omega(6)];

if (omega(1) < 0 || omega(3) < 0 || omega(6) < 0)
    w = -w;
end

% k = chol(w);
% K = ((inv(k))'*k(3,3))';

%intrinsic parameters 
y0 = (omega(2)*omega(4)-omega(1)*omega(5))/(omega(1)*omega(3)-omega(2)^2);
lambda = omega(6)-(omega(4)^2+y0*(omega(2)*omega(4)-omega(1)*omega(5)))/omega(1);
ax = sqrt(lambda/omega(1));
ay = sqrt(lambda*omega(1)/(omega(1)*omega(3)-omega(2)^2));
s = -omega(2)*ax^2*ay/lambda;
x0 = s*y0/ay-omega(4)*ax^2/lambda;
K = [ax s x0; 0 ay y0; 0 0 1];

%extrinsic parameters
p = zeros(1,11);
p(1:5) = [ax s x0 ay y0];
Kinv = inv(K);
t = Kinv*H(:,3);
mag = norm(Kinv*H(:,1));
if(t(3)<0)
mag = -mag;
end
r1 = Kinv*H(:,1)/mag;
r2 = Kinv*H(:,2)/mag;
r3 = cross(r1,r2);
R = [r1 r2 r3];
t = t/mag;
[U,D,V] = svd(R);
R = U*V';

P = K*[R', t];

% camera coordinates 
Xc = -det([P(:, 2), P(:, 3), P(:, 4)]);
Yc = det([P(:, 1), P(:, 3), P(:, 4)]);
Zc = -det([P(:, 1), P(:, 2), P(:, 4)]);
Wc = det([P(:, 1), P(:, 2), P(:, 3)]);
CC = [Xc, Yc, Zc, Wc];
CC = CC/CC(4);

% camera settings 
%orientation = R';
%location = -t' * orientation;

% first we find the scale for the measure we know 
AB = 243;
a_y = ((P(:, 2)/P(3, 2) - b2) \ (b2 - a2))/AB;

% discover the distance of the line in common 
horizon = real(cross(P(:, 3)/P(3, 3), P(:, 1)/P(3, 1)));
length = sqrt(horizon(1)^2 + horizon(2)^2);
horizon = horizon/length;

line1 = real(cross(a2, d2));
v = real(cross(line1, horizon));
v = v/v(3); 

line2 = real(cross(v', b2));
vertical_line = real(cross(a2, d2)); 
t = real(cross(line2, vertical_line));
t = t/t(3);
AD = AB*sqrt(sum((a2-d2).^2))*sqrt(sum(((P(:,2)/P(3, 2))'-t).^2))/...
    sqrt(sum((t'-d2).^2))/sqrt(sum(((P(:,2)/P(3, 2))-a2).^2));


% same procedure 
a_x = ((P(:, 1)/P(3, 1) - d2) \ (d2 - a2))/AD;

horizon = real(cross(P(:, 3)/P(3, 3), P(:, 2)/P(3, 2)));
length = sqrt(horizon(1)^2 + horizon(2)^2);
horizon = horizon/length;

line1 = real(cross(d2, h2));
v = real(cross(line1, horizon));
v = v/v(3); 

line2 = real(cross(v', a2));
vertical_line = real(cross(d2, h2)); 
t = real(cross(line2, vertical_line));
t = t/t(3);
DH = AB*sqrt(sum((d2-h2).^2))*sqrt(sum(((P(:,1)/P(3, 1))'-t).^2))/...
    sqrt(sum((t'-h2).^2))/sqrt(sum(((P(:,1)/P(3, 1))-d2).^2));

% and finally 
a_z = ((P(:, 3)/P(3, 3) - h2) \ (h2 - d2))/DH;

P = [P(:,1)*a_x, P(:,2)*a_y, P(:,3)*a_z, P(:,4)];
P = [P(:,1) -P(:,2) -P(:,3) P(:,4)];

Hxy = [P(:,1:2),P(:,4)];
Hxz = [P(:,1),P(:,3:4)];
Hyz = P(:,2:4);

Halt = [P(:,1), P(:,2), P(:,4)];
worldPoints(1, :) = Halt * a;
worldPoints(2, :) = Halt * b;
worldPoints(3, :) = Halt * c;
worldPoints(4, :) = Halt * d;
worldPoints(5, :) = Halt * a2;
worldPoints(6, :) = Halt * b2;
worldPoints(7, :) = Halt * c2;
worldPoints(8, :) = Halt * d2;
% worldPoints(9, :) = Halt * e2;
% worldPoints(10, :) = Halt * h2;
% worldPoints(11, :) = Halt * e;
% worldPoints(12, :) = Halt * h;
worldPoints(:, 4) = 1;
