mytitle = 'solid body rotation test';
lag_pos_file = 'solidbody.dat';

x_release = 250.05:100:750.05; nlag = numel(x_release);
y_release = 500*ones(nlag,1);
t_release = greg2mjulian(2009,1,1,0,0,0);  %release on april 1, 2007
s_release = -rand(nlag,1);                 %release at various heights
gid       = zeros(nlag,1);

%--------------------------------------------------------------
% dump to ascii file 
%--------------------------------------------------------------
write_offline_lagfile_ascii(lag_pos_file,mytitle,x_release, ...
    y_release,s_release,'s',t_release*ones(nlag,1));

% -------------------------------------------------------
% need a post-processing ascii file
% want daily position of particles and averaged daily depth
%----------------------------------------------------------
