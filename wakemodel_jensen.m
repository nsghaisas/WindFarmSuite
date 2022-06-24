%clear
function [uq, u_diskavg, u_turbavg_zh, u_horzavg_zh] =  wakemodel_jensen(ntx, nty, nturb, xturb, yturb, zturb, zh, diam, ct, kwinp, ustarinf, z0lo, kap, iquery_points, compute_full_x, xg, yg, zg, outdir)

  %%---------beg inputs block----------------------------------------------------------------------
  %ntx = 12; nty = 6;
  %Lx = 3*pi; Ly = pi; Lz = 1;
  %zh = 0.1; diam = 0.1;
  %
  %sx = Lx/ntx; sy = Ly/nty;
  %
  %ustarhi = 1; z0lo = 1e-4;
  %ct = 3/4; kap = 0.4; kwinp = 0.035;
  %
  %nturb = ntx*nty;
  %xturb  = zeros(nturb,1); yturb  = zeros(nturb,1); zturb  = zeros(nturb,1);
  %for itx = 1:ntx
  %  for jty = 1:nty
  %    it = (itx-1)*nty + jty;
  %    xturb(it) = (itx-0.5)*sx;
  %    yturb(it) = (jty-0.5)*sy;
  %  end
  %end
  %zturb(:) = zh;
  %
  %iquery_points = 1;
  %%---------end inputs block----------------------------------------------------------------------
  
  betturb = zeros(nturb,1); ctturb = zeros(nturb,1); kwturb = zeros(nturb,1);
  ctturb(:) = ct;
  betturb = (1 + sqrt(1-ctturb))./(2*sqrt(1-ctturb));
  %betturb(:) = 1;
  kwturb(:) = kwinp;
  
  irow_1st = 1:nty:nturb;    icol_1st = 1:1:nty;
  irow_mid = ceil(nty/2):nty:nturb;
  
  % query points
  if(iquery_points==1)
      % query at turbine locations
  
      ndy = 100; yline = linspace(0, diam, ndy);
      ndz = 100; zline = linspace(0, diam, ndz);
      [y2d, z2d] = ndgrid(yline, zline);
  
      xq = zeros(nturb, ndy, ndz);
      yq = zeros(nturb, ndy, ndz);
      zq = zeros(nturb, ndy, ndz);
      for it = 1:nturb
        xq(it,:,:) = xturb(it);
        yq(it,:,:) = y2d + (yturb(it)-diam/2);
        zq(it,:,:) = z2d + (zturb(it)-diam/2);
      end

  elseif(iquery_points==2)
      % query on xg, yg, zg
      xq = xg; yq = yg; zq = zg;
      uq = zeros(size(xq));
  end
  
  mindist = 0.001*diam;
  maskloc = zeros(size(xq(1,:,:)));
  maskdsk = zeros(size(xq(1,:,:)));
  sumloc  = zeros(size(xq(1,:,:)));
 
  % compute area_frac for averaging over disk
  if(iquery_points==1)
      r2d = sqrt((y2d-diam/2).^2 + (z2d-diam/2).^2);
      maskdsk(r2d < diam/2) = 1;
      area_frac = sum(maskdsk,'all')/length(maskdsk(:));
      %[sum(maskdsk,'all')   length(maskdsk(:))    area_frac]
  end
 
  uinf = zeros(size(xq(1,:,:)));  %zeros(1,ndy, ndz);
  %size(uinf)
  %size(xq)
  %size(xq(1,:,:))
  for j=1:size(xq,2)
    uinf(1,j,:) = ustarinf/kap*log(zq(1,j,:)/z0lo);
  end

  if(iquery_points==1)
      if(compute_full_x==1) 
        i_arr = irow_mid;
      else
        i_arr = [irow_mid(end-2) irow_mid(end-1) irow_mid(end)];
      end
  elseif(iquery_points==2)
      if(compute_full_x==1) 
        ixqst = 1; 
      else
        frac_compute_inx = ceil(length(irow_mid)/3);  % last 3 rows of the middle column
        %frac_compute_inx = 10;
        ixqst = ceil(size(xq,1)*(1-1/frac_compute_inx)); %ixqst
      end
      i_arr = ixqst:size(xq,1);
  end
  for i = i_arr
  %if(compute_full_x==1) 
  %  ixqst = 1; 
  %else
  %  frac_compute_inx = ceil(length(irow_mid)/3);  % last 3 rows of the middle column
  %  %frac_compute_inx = 10;
  %  ixqst = ceil(size(xq,1)*(1-1/frac_compute_inx)); %ixqst
  %end
  %ixqst
  %size(xq,1)
  %for i = ixqst:size(xq,1)
    sumloc(:) = 0;
    for iup = 1:nturb
      xdist = xq(i) - xturb(iup);
  
      if(xdist >  mindist)
          maskloc(:) = 0;
          raddist = sqrt((yq(i,:,:)-yturb(iup)).^2 + (zq(i,:,:)-zturb(iup)).^2);
          maskloc(raddist < (diam/2 + kwturb(iup)*xdist)) = 1;
          du_one = betturb(iup)*(1-sqrt(1-ctturb(iup)))/(sqrt(betturb(iup)) + kwturb(iup)*xdist/(diam/2))^2;
          sumloc = sumloc + du_one^2*maskloc;
      end
    end
    uq(i,:,:) = ones(size(sumloc)) - sqrt(sumloc);
    uq(i,:,:) = uq(i,:,:) .* uinf(1,:,:);
  
  end

  % determine interpolation factors in vertical
  [vv, izh] = min(abs(zq(1,1,:)-zh));
  if(zq(1,1,izh) > zh)
    i1 = izh-1; i2 = izh;
  else
    i1 = izh; i2 = izh+1;
  end

  %size(zq)
  %[i1 i2]
  %zh
  %irow_mid

  u_diskavg    = zeros(length(irow_mid),1); 
  u_turbavg_zh = zeros(length(irow_mid),1); 
  u_horzavg_zh = zeros(length(irow_mid),1);
  if(iquery_points==1)
      udisk = mean(uq.*maskdsk, [2 3])/area_frac;
      %(udisk(irow_1st)'./udisk(1)).^3
      u_diskavg(:) = udisk(irow_mid);

      icounter = 0;
      for it = irow_mid
          uplane = uq(it,:,i1) + (zh-zq(1,1,i1))/(zq(1,1,i2)-zq(1,1,i1)) * (uq(it,:,i2) - uq(it,:,i1));
          icounter = icounter+1;
          u_turbavg_zh(icounter) = mean(uplane);
          %[it i1 i2 zh zq(1,1,i1) zq(1,1,i2) size(uplane)]

          %figure(1), clf
          %contourf(squeeze(yq(it,:,:)), squeeze(zq(it,:,:)), squeeze(uq(it,:,:)),'edgecolor','none'), daspect([1 1 1])
          %xlabel('y'), ylabel('z'), set(gca,'fontsize',14)
          %colormap jet
          %caxis([10 18])
          %colorbar
          %fname = strcat(outdir,'/ucontours_iq01_',sprintf('%4.4i',size(xg,1)),'_turb',sprintf('%3.3i',it),'.png');
          %screen2jpeg(fname)
      end
  elseif(iquery_points==2)
      sx = xturb(nty+1)-xturb(1);
      icounter = 0;
      for it = irow_mid
        icounter = icounter+1;
        %[it irow_mid(end-1) irow_mid(end)]
        if((compute_full_x==1) || (it==irow_mid(end)) || (it==irow_mid(end-1)))
          xst = xturb(it) - sx/2; xen = xturb(it) + sx/2;
          [vv,ist] = min(abs(xq(:,1,1)-xst));   [vv,ien] = min(abs(xq(:,1,1)-xen));
          uplane = uq(ist:ien,:,i1) + (zh-zq(1,1,i1))/(zq(1,1,i2)-zq(1,1,i1)) * (uq(ist:ien,:,i2) - uq(ist:ien,:,i1));
          u_horzavg_zh(icounter) = mean(uplane, [1 2]);
        end
      end
      %figure(1), clf
      %contourf(squeeze(xq(:,1,1)), squeeze(yq(1,:,1)), squeeze(uq(:,:,izh))','edgecolor','none'), daspect([1 1 1])
      %xlabel('x'), ylabel('y'), set(gca,'fontsize',14)
      %colormap jet
      %%caxis([0.8 2.2])
      %colorbar
      %grid on
      %fname = strcat(outdir,'/ucontours_iq02_',sprintf('%4.4i',size(xq,1)),'.png');
      %screen2jpeg(fname)
  end
  %xq(1:5,1:3,1)
  %yq(1:5,1:3,1)

  %% check if farm is large enough
  %u_horzavg_diff = 1.0  - u_horzavg_zh(end-1) / u_horzavg_zh(end);
  %if(abs(u_horzavg_diff)>0.01)
  %   'Farm is not large enough: '
  %   u_horzavg_diff
  %end

  %% check if farm is large enough
  %u_turbavg_diff = 1.0  - u_turbavg_zh(end-1) / u_turbavg_zh(end);
  %if(abs(u_turbavg_diff)>0.01)
  %   'Farm is not large enough: '
  %   u_turbavg_diff
  %end
end
