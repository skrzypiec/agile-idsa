  function agile_idsa(cyc,variable,range,col)

% function agile_idsa(cyc,variable,range,col)
%
% cyc ... cycle number for first run
% variable = 1 (radius) or 2 (rest mass)
% range ... range of independent variable to display
% col ... color

% read agile data
  [status,result] = unix('ls -1 *.stp|head -n1');
  len = size(result,2);
  name = result(1:len-9);
  zstr = '0000';
  cycstr = int2str(cyc);
  cycstr = [zstr(1:4-size(cycstr,2)),cycstr];
  filename = [name,cycstr,'.stp'];
  disp(['reading ',filename]);
  [nq,ny,t,label,y,dx,s,cs,status] = getagile(filename);

% read nuprox data
  zstr = '00000';
  cycstr = int2str(cyc);
  cycstr = [zstr(1:5-size(cycstr,2)),cycstr];
  filename = ['nuprox1d',cycstr,'.d'];
  disp(['reading ',filename]);
  [nr,nd,nt,nye,nmue,nmuhat,nynu,nznu,ntnu,neta,...
    nyedot,nedot,nLN,nLE] = getnuprox(filename,nq-1);

%.prepare data..........................................................
  Ms = 1.9890e+33;
  amr = y(:,1)/Ms;
  ar = y(:,2);
  av = y(:,3);
  bd = y(2:nq,5);
  bye = y(2:nq,7);
  bs = s(2:nq);
  nLN(2:nq,:) = nLN;
  nLN(1,:) = 0;
  nLE(2:nq,:) = nLE;
  nLE(1,:) = 0;

%.independent variable
  if variable==1
    ax = ar;              %left edge
    bx(1:nq-1) = (0.5*(ar(1:nq-1).^3+ar(2:nq).^3)).^(1/3);
    xlab = 'Radius [cm]';
  else
    ax = amr;             %left edge
    bx(1) = amr(1);
    bx(1:nq-1) = 0.5*(amr(1:nq-1)+amr(2:nq));
    xlab = 'Mass [Ms]';
  end

% velocities
  subplot(3,3,1);
  hold on;
  box on;
  h = plot(ax,av,'Color',col);
  axis('auto');
  v = axis;
  v(1:2) = range;
  axis(v);
  xlabel(xlab);
  ylabel('Velocity [cm/s]');

%.lepton fraction
  subplot(3,3,5);
  hold on;
  box on;
  nyl = nye + nynu(:,1) - nynu(:,2);
  h = plot(bx,nyl,'Color',col);
  axis('auto');
  v = axis;
  v(1:2) = range;
  axis(v);
  xlabel(xlab);
  ylabel('YL []');

%.lepton number luminosity
  subplot(3,3,4);
  hold on;
  box on;
  h1 = plot(ax,nLN(:,1),'LineStyle','-','Color',col);
  h2 = plot(ax,nLN(:,2),'LineStyle','--','Color',col);
  axis('auto');
  v = axis;
  v(1:2) = range;
  axis(v);
  xlabel(xlab);
  ylabel('Lepton Lum. [#/s]');

%.density
  subplot(3,3,3);
  hold on;
  box on;
  h = plot(bx,bd,'Color',col);
  axis('auto');
  v = axis;
  v(1:2) = range;
  axis(v);
  xlabel(xlab);
  ylabel('Density [g/cm^3]');
  set(gca,'YScale','log','YTick',[1e6 1e8 1e10 1e12 1e14]);

%.electron fraction
  subplot(3,3,2);
  hold on;
  box on;
  h = plot(bx,bye,'Color',col);
  axis('auto');
  v = axis;
  v(1:2) = range;
  axis(v);
  xlabel(xlab);
  ylabel('Ye []');

%.entropy
  subplot(3,3,6);
  hold on;
  box on;
  h = plot(bx,bs,'Color',col);
  ylabel('Entropy [kB/baryon]');
  axis('auto');
  v = axis;
  v(1:2) = range;
  axis(v);
  xlabel(xlab);

%.neutrino luminosity
  subplot(3,3,7);
  hold on;
  box on;
  h1 = plot(ax,nLE(:,1),'Color',col,'LineStyle','-');
  h2 = plot(ax,nLE(:,2),'Color',col,'LineStyle','--');
  axis('auto');
  v = axis;
  v(1:2) = range;
  axis(v);
  xlabel(xlab);
  ylabel('Energy Lum. [erg/s]');

%.Ynu
  subplot(3,3,8);
  hold on;
  box on;
  h1 = plot(bx,nynu(:,1),'Color',col,'LineStyle','-');
  h2 = plot(bx,nynu(:,2),'Color',col,'LineStyle','--');
  axis('auto');
   v = axis;
   v(1:2) = range;
   axis(v);
  xlabel(xlab);
  ylabel('Y\nu []');

%.Znu
  subplot(3,3,9);
  hold on;
  box on;
  nynu = nynu + 1.d-08;
  h1 = plot(bx,nznu(:,1)./nynu(:,1),'Color',col,'LineStyle','-');
  h2 = plot(bx,nznu(:,2)./nynu(:,2),'Color',col,'LineStyle','--');
  axis([range,0,50]);
  xlabel(xlab);
  ylabel('Energy/neutrino [MeV]');
