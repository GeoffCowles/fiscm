z =0:.5:40;
kh = .01;

fid = fopen('kh_data3.dat','w');
for i=1:numel(z);
  fprintf(fid,'%f %f\n',z(i),kh);
end;
fclose(fid);
