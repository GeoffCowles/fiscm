
clear all
%close all
fname = 'fiscm_1.nc';
delay = .2;
figure


% open the particle data
nc = netcdf(fname,'nowrite');

% read particl position data
zp = nc{'z'}(:,:);
hp = nc{'h'}(:,:);
status = nc{'status'}(:,:);


% determine problem size dimensions
dims = size(zp);
ntimes = dims(1); 
nlag   = dims(2);



% initial locations
%plot(xp(1,:),yp(1,:),'ro','MarkerFaceColor','r','EraseMode','none');
%plothandle = plot(xp(1),yp(1),'bo','MarkerFaceColor','b','EraseMode','xor');
%figure(fighandle)

for n=1:1:ntimes
  for i=nlag:nlag %1:1 %nlag
    if(status(n,i)>0); plot(n,zp(n,i),'k+');  hold on; end;
    if(status(n,i)<0); plot(n,zp(n,i),'r+');  hold on; end;
    plot(n,-hp(n,i),'b+');  hold on; 
  end
  axis([1,ntimes,-100,2])
% set(plothandle,'XData',xp(n,:),'YData',yp(n,:))
% drawnow
 pause(delay)
end;

