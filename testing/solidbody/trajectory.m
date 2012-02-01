clear all
close all
fname = '../../fiscm_1.nc';
delay = .1;
stride=1;

% open the mesh file
nc = netcdf('solidbody.nc','nowrite');
xm = nc{'x'}(:);
ym = nc{'y'}(:);
hm = nc{'h'}(:);
nv = nc{'nv'}(:)';

fighandle = figure;
patch('Vertices',[xm,ym],'Faces',nv,...
       'Cdata',hm,'edgecolor','interp','facecolor','interp');
axis([0 1000 0 1000]);
hold on

% open the particle data
nc = netcdf(fname,'nowrite');

% read particl position data
time = nc{'time'}(:);
xp = nc{'x'}(:,:);
yp = nc{'y'}(:,:);
%up = nc{'u'}(:,:);
%vp = nc{'v'}(:,:);
%tp = nc{'T'}(:,:);
%cp = nc{'cell'}(:,:);


i = 1; xplot=xp(i,:); yplot=yp(i,:);
h = plot(xplot,yplot,'k.','EraseMode','none');
for i=1+stride:stride:numel(time)
    drawnow; 
    xplot=xp(i,:); yplot=yp(i,:);
    plot([xp(i,:);xp(i-stride,:)],[yp(i,:);yp(i-stride,:)],'r-')
    set(h,'xdata',xplot,'ydata',yplot)
    pause(delay)

end
plot(xp(end,:),yp(end,:),'b.','EraseMode','none');