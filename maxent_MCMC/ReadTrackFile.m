function track = ReadTrackFile(dataFile)
%keyboard
fid = fopen(dataFile, 'r');
track.converged=false;
try
    % retrieve track log
    A = load(dataFile);
    track.all=A;
    track.outer=int16(A(:,1));
    track.inner=int16(A(:,2));
    track.metrics=A(:,3:5);
    a=diff(track.metrics);
    a=a(a(:,1)~=0,:);
    track.metricsteps=a;
    track.total_time=A(:,6);
    track.sample_time=A(:,7);
    track.learn_time=A(:,8);
    track.loss_ent=A(:,9:10);
    track.exit_conds=logical(A(:,11:14));
    fclose(fid);

    if size(track.exit_conds,1)<=1
        track.exit='no trackfile';
    elseif all(track.exit_conds(end-1,1:3)==1)
        track.exit='reached target';
        track.converged=1;
    elseif track.exit_conds(end-1,end)==1;
        track.exit='time ran out';
    else
        track.exit='maximal number of runs reached, or something strange happened..';
    end

    track.header = ['Outer Loop | Inner Loop | Metric mean | Metric corr | Metric cov | time | stat1 | stat2 | exit_conditions '];
    track.success=true;
catch
    track.exit='track-file not found';
    track.success=false;
end