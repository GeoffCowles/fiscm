function write_offline_lagfile_ascii(file,info,x,y,v,vtype,trel)
% write position and release time data for particles to initialize offline Lagrangian code
%
% function write_offline_lagfile_ascii(info,info,x,y,s,trel)
%
% DESCRIPTION:
%   write position and release time data for particles to initialize 
%   offline Lagrangian code
%
% INPUT
%   file:  filename
%   info:  string containing information about case (a global attribute)
%   x:     x-coordinate of initial particle position
%   y:     y-coordinate of initial particle position
%   v:     vertical-coordinate of initial particle position
%   vtype: vertical coordinate type (='z' for z-coordinate, ='s' for sigma coordinate)
%
% OUTPUT:
%    Initial particle position file
%
% Author(s):  
%    Geoff Cowles (University of Massachusetts Dartmouth)
%
% Revision history
%   
%==============================================================================

fprintf('creating initial lagrangian position file\n')

% open boundary forcing
fid = fopen(file,'w');
fprintf(fid,'%d %d\n',numel(x),2);
fprintf(fid,'%s\n',info);
if(vtype == 'z')
  fprintf(fid,'vertical coordinate = z');
  fprintf(fid,'x     y     z    t_release\n');
else
  fprintf(fid,'x     y     s    t_release\n');
end;

for i=1:numel(x)
  fprintf(fid,'%10.1f \t %10.1f \t %f \t %f\n',x(i),y(i),v(i),trel(i));
end;
fclose(fid);

