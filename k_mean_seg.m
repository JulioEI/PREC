function [Ik, k] = k_mean_seg(I1, varargin)


sj = inputParser;
addParameter(sj,'nk',5,@isnumeric)
addParameter(sj,'v_range',[0 240],@isnumeric)
parse(sj,varargin{:})
nk = sj.Results.nk;
v_range = sj.Results.v_range;

    [r, c] = size(I1);
    I1 = reshape(I1,1,r*c);
    I1 = double(I1);
    Ik = double(zeros(size(I1)));
    if isempty(v_range)
        v_range = [0 240];
    end
    k = linspace(v_range(1), v_range(2), nk);
    kt = zeros(2,nk);
    
    it = 0;
    stop = 0;
    max_it = 20;
    while((stop == 0)&&(it<max_it))
        stop = 1;
        kt(:,:) = 0;
        %EXPECTATION - Assign cluster to each point 
        Ikt = double(I1);
        for i = 2:nk
            Ikt= double([Ikt; I1]);
        end
        Ikt = Ikt - k';
        [~, cluster] = min(sqrt(Ikt.^2));
        for i = 1:size(cluster,2)
            value = cluster(i);
            kt(1,value) = kt(1,value) + I1(1,i);
            kt(2,value) = kt(2,value) + 1;
            Ik(1, i) = k(value);
        end
        %MAXIMIZATION - Move clusters centroid to the mean of the assigned points
        a = kt(1,:)./(kt(2,:));
        if (abs(max(abs(a-k)))>0.1)
            stop=0;
            k=a;
        end
        it = it+1;
    end
    Ik = reshape(Ik, r ,c);
end



