clear all
close all


fname = 'fake_forcing.nc';
ncid = netcdf.open(fname);
kh = squeeze(netcdf.getVar(ncid,12,[1,0,1],[1,81,1]));
nlay = numel(kh);
depth = 40;
z = 0:-(40-0)/(nlay-1):-40;
subplot(1,2,1)
plot(kh,z);
xlabel('kh');
ylabel('z');


fname = 'fiscm_1.nc';
delay = .2;


% open the particle data
subplot(1,2,2)
ncid = netcdf.open(fname);
[numdims,numvars,numglobalatts,unlimdimid] = netcdf.inq(ncid);

% read particle position data
time = netcdf.getVar(ncid,0);
zp = netcdf.getVar(ncid,7);
hp = netcdf.getVar(ncid,3);


% determine problem size dimensions
[nlag,ntimes] = size(zp);


% initial locations
%plot(xp(1,:),yp(1,:),'ro','MarkerFaceColor','r','EraseMode','none');
%plothandle = plot(xp(1),yp(1),'bo','MarkerFaceColor','b','EraseMode','xor');
%figure(fighandle)
prob = zeros(40,ntimes); 
for n=1:ntimes
  [nbin] = histc(abs(zp(:,n)),0:1:40);
%  tt = ceil(abs(zp(:,n))+1e-9);     
%  [nbin,x] = hist(abs(zp(:,n)),40);
  prob(40:-1:1,n) = 40*nbin(1:end-1)/nlag; 
%  plot(n*ones(nlag,1),zp(:,n),'k+'); hold on;
%  axis([1,ntimes,-100,2])
% set(plothandle,'XData',xp(n,:),'YData',yp(n,:))
% drawnow
% pause(delay)
end;
pcolor(prob); shading interp; colorbar; caxis([0,2])

