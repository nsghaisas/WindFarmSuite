clear

% set wind farm coordinates
% xt, yt, zt, diam, CT
turbine_data =   [1000 1000 100 100 0.75;
                  1700 1000 100 100 0.75;
                  2400 1000 100 100 0.75;
                  3100 1000 100 100 0.75;
                  3800 1000 100 100 0.75;
                  1000 1400 100 100 0.75;
                  1700 1400 100 100 0.75;
                  2400 1400 100 100 0.75;
                  3100 1400 100 100 0.75;
                  3800 1400 100 100 0.75;
                  1000 1800 100 100 0.75;
                  1700 1800 100 100 0.75;
                  2400 1800 100 100 0.75;
                  3100 1800 100 100 0.75;
                  3800 1800 100 100 0.75;
                  1000 2200 100 100 0.75;
                  1700 2200 100 100 0.75;
                  2400 2200 100 100 0.75;
                  3100 2200 100 100 0.75;
                  3800 2200 100 100 0.75];


ntx = 5; nty = 4; nturb = ntx*nty;
xt = turbine_data(:,1); yt = turbine_data(:,2); zt = turbine_data(:,3);
DIAM = turbine_data(:,4); CT = turbine_data(:,5);

% set other inputs; 
diam = DIAM(1); ct = CT(1);  % since diam and ct are same for each turbine
kwinp = 0.035;               % wake expansion coefficient
zplane = 100;                % what plane do you want the results?
ustarinf = 0.45;             % define upstream profile friction velocity
z0lo = 0.1;                  % surface roughness of the ground
kap = 0.41;                  % von-Karman constant for turbulent boundary layers
iquery_points = 2;           % 2 :: compute on a grid; 1 :: compute on turbine disks
compute_full_x = 1;          % 1 :: compute everywhere; 2 :: only last 3 rows to speed up calculation
dx=10; dy=10; dz=50;         % resolution in x,y,z
winddir = 315;               % what direction is the wind coming from
outdir = './outputs/';       % where do you want to write outputs
mkdir(outdir)

% rotate turbines according to wind direction
thet_rot = winddir-270;
xref = 0.5*(min(xt)+max(xt));
yref = 0.5*(min(yt)+max(yt));
for ii=1:nturb
  xl = xt(ii)-xref; yl = yt(ii)-yref;
  xrot = cos(thet_rot*pi/180)*xl - sin(thet_rot*pi/180)*yl;
  yrot = sin(thet_rot*pi/180)*xl + cos(thet_rot*pi/180)*yl;
  xt(ii) = xrot + xref; yt(ii) = yrot + yref;
end

% set domain extents
xl = min(xt) - 15*max(diam);    xr = max(xt) + 15*max(diam);
yl = min(yt) - 15*max(diam);    yr = max(yt) + 15*max(diam);
zl = 0.0;                       zr = 1000;

% construct the grid
xline = xl:dx:xr; yline = yl:dy:yr; zline = zl:dz:zr;
[xg, yg, zg] = ndgrid(xline, yline, zline);

% evaluate wake model
[uq, u_diskavg, u_turbavg_zh, u_horzavg_zh] =  wakemodel_jensen(ntx, nty, nturb, xt, yt, zt, zplane, diam, ct, kwinp, ustarinf, z0lo, kap, iquery_points, compute_full_x, xg, yg, zg, outdir);

% plot contours at the desired plane
[vv, izh] = min(abs(zg(1,1,:)-zplane));
figure(1), clf
contourf(squeeze(xg(:,1,1)), squeeze(yg(1,:,1)), squeeze(uq(:,:,izh))','edgecolor','none'), daspect([1 1 1])
xlabel('x'), ylabel('y'), set(gca,'fontsize',14)
colormap jet, colorbar
hold on
for ii=1:nturb
  plot(xt(ii),yt(ii),'.k','MarkerSize',20)
end
fname = strcat(outdir,'/ucontours_iq',sprintf('%2.2i', iquery_points),'_dir',sprintf('%3.3i',winddir),'_',sprintf('%4.4i',size(xg,1)),'.png');
screen2jpeg(fname)

