
clear all
close all
fname = 'fiscm_1.nc';
delay = .2;

% open the mesh file
nc = netcdf('mesh.nc','nowrite');
xm = nc{'x'}(:);
ym = nc{'y'}(:);
hm = nc{'h'}(:);
nv = nc{'nv'}(:)';


fighandle = figure;
coords(:,1) = xm;
coords(:,2) = ym;
patch('Vertices',coords,'Faces',nv,...
       'Cdata',hm,'edgecolor','interp','facecolor','interp');
colorbar
hold on

% open the particle data
nc = netcdf(fname,'nowrite');
times = nc{'time'}(:);

% read particl position data
time = nc{'time'}(:);
xp = nc{'x'}(:,:);
yp = nc{'y'}(:,:);
up = nc{'u'}(:,:);
vp = nc{'v'}(:,:);
tp = nc{'T'}(:,:);
cp = nc{'cell'}(:,:);

%figure
%pt = 20;
%velmag = sqrt(  up(:,pt).^2 + vp(:,pt).^2)   ;
%plot(time,velmag)
%figure
%plot(xp(:,pt),yp(:,pt))
%%pause



% determine problem size dimensions
dims = size(xp);
ntimes = dims(1); 
nlag   = dims(2);

% compute the plot limits
xmax = max(xp(:));
xmin = min(xp(:));
ymax = max(yp(:));
ymin = min(yp(:));
dx = xmax-xmin;
dy = ymax-ymin
fac = .1;
%axis([xmin-fac*dx,xmax+fac*dx,ymin-fac*dy,ymax+fac*dy])
hold on;
%axis([1.1e6,1.3e6,-2e5,0.5e5])


% initial locations
plot(xp(1,:),yp(1,:),'ro','MarkerFaceColor','r','EraseMode','none');
plothandle = plot(xp(1),yp(1),'bo','MarkerFaceColor','b','EraseMode','xor');
figure(fighandle)

for n=2:ntimes %ntimes
  for i=1:nlag %nlag 
    xseg = [xp(n,i),xp(n-1,i)];
    yseg = [yp(n,i),yp(n-1,i)];
%   plot(xseg,yseg,'k-','EraseMode','none'); hold on;
%ok   plot(xseg,yseg,'w-'); hold on;
   if(i==1)
   plot(xseg,yseg,'k-'); hold on;
   elseif(i==2)
   plot(xseg,yseg,'r-'); hold on;
   else
   plot(xseg,yseg,'w-'); hold on;
   end;  
  end;
% set(plothandle,'XData',xp(n,:),'YData',yp(n,:))
 drawnow
 n
 pause(delay)
end;


