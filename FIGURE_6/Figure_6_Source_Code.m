%
%      MATLAB Script for Fig. 6 in the eLife manuscript:
%
%  "Nanophysiology Approach Reveals Diversity in Ca2+ Microdomains ..."
%   Rameshkumar, Shrestha, Boff, Hoon, Matveev, Zenisek, Vaithianathan
%   https://elifesciences.org/reviewed-preprints/105875#s2
%               Code: Victor Matveev, Sep 1, 2025
% ==============================================================
% 3D Geometry Visualization: Box, Ribbon, Stalk, Section Plane
% ==============================================================

% ----- Geometry parameters -----------------------------------
zAxis    = 0.190;
xAxis    = zAxis;
yAxis    = zAxis * 0.375;
rStalk   = 0.015;
zElevate = 0.03;
z0       = zAxis + zElevate;

Lz = 1.10;  Lx = 0.64;  Ly = Lx;

% ----- Figure setup ------------------------------------------
figure; hold on; view(-37, 30);

% ==============================================================
%                       Draw Box Frame
% ==============================================================
xx = [-Lx  Lx  Lx -Lx -Lx];
yy = [-Ly -Ly  Ly  Ly -Ly];
plot3(xx, yy,  0*[1 1 1 1 1], 'k-', 'LineWidth', 2);
plot3(xx, yy, Lz*[1 1 1 1 1], 'k-', 'LineWidth', 2);
plot3([-Lx -Lx], [-Ly -Ly], [0 Lz], 'k-', 'LineWidth', 2);
plot3([-Lx -Lx], [ Ly  Ly], [0 Lz], 'k-', 'LineWidth', 2);
plot3([ Lx  Lx], [-Ly -Ly], [0 Lz], 'k-', 'LineWidth', 2);
plot3([ Lx  Lx], [ Ly  Ly], [0 Lz], 'k-', 'LineWidth', 2);

% ==============================================================
%                Draw Ribbon (capsule surfaces)
% ==============================================================
N = 40;
rr   = linspace(0,1,N);
thth = linspace(0,2*pi,N);
[r, th] = meshgrid(rr, thth);

x = xAxis*r.*cos(th);
y = yAxis*r.*sin(th);
zTop = z0 + zAxis*sqrt(1 - r.^2);
zBot = z0 - zAxis*sqrt(1 - r.^2);

surfl(x, y, zTop, [-45 30], [.65 .4 .3 10]);
surfl(x, y, zBot, [135 30], [.65 .4 .3 10]);
colormap(bone); shading flat;

% ==============================================================
%                Draw Stalk (arciform density)
% ==============================================================
xx = 0.03; yy = 0.015; zz = 0.05;
coord = [...
   -xx -yy  0;  xx -yy  0;  xx  yy  0; -xx  yy  0; ...
   -xx -yy zz;  xx -yy zz;  xx  yy zz; -xx  yy zz];
idx = [4 8 5 1 4; 1 5 6 2 1; 2 6 7 3 2; ...
       3 7 8 4 3; 5 8 7 6 5; 1 4 3 2 1]';
patch(coord(idx,1), coord(idx,2), coord(idx,3), ...
      [0.508 0.577 0.633], 'FaceAlpha', 1);

% ==============================================================
%                 Draw Section Plane (x = 0)
% ==============================================================
[y, z] = meshgrid([-Ly Ly], [0 Lz]);
x = zeros(size(y));
surf(x, y, z, ...
     'FaceAlpha', 0.3, ...            % semi-transparent
     'FaceColor', [0.4 0.4 0.6], ...  % darker bluish tint
     'EdgeColor', 'none');            % cleaner look

% ==============================================================
%                    Draw Channel Disks
% ==============================================================
ICa_X = rStalk + 0.025;
xx = [-ICa_X, ICa_X, -ICa_X, ICa_X];
yy = [-ICa_X, -ICa_X, ICa_X, ICa_X];
[rr, th] = meshgrid(linspace(0.01,0.02,10), linspace(0,2*pi,30));

for k = 1:numel(xx)
    x0 = xx(k); y0 = yy(k);
    surf(x0 + rr.*cos(th), y0 + rr.*sin(th), zeros(size(rr)));
end

% ==============================================================
%                        Axes, Labels
% ==============================================================
set(gca, 'XTick', [-Lx 0 Lx], ...
         'YTick', [-Ly 0 Ly], ...
         'ZTick', [0 Lz], ...
         'FontSize', 14);
xlabel('Y, \mum', 'FontSize', 11);
ylabel('X, \mum', 'FontSize', 11);
zlabel('Z, \mum', 'FontSize', 11);
axis tight;

