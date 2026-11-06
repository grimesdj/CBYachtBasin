function [Pnear, dist, iseg, t, ivert] = nearest_point_on_curve(P, V)
%NEARESTPOINTONCURVE Find closest point on a piecewise-linear curve.
%
% INPUTS
%   P    - N x 2 array of query points [x y]
%   V    - M x 2 array of curve vertices [x y]
%
% OUTPUTS
%   Pnear - N x 2 nearest points on the curve
%   dist  - N x 1 Euclidean distance to the curve
%   iseg  - N x 1 index of the nearest curve segment
%          (segment is V(iseg,:) -> V(iseg+1,:))
%   t     - N x 1 normalized position along nearest segment
%          t = 0 at V(iseg,:)
%          t = 1 at V(iseg+1,:)
%   ivert - N x 1 index of nearest vertex

    % Check dimensions
    if size(P,2) ~= 2 || size(V,2) ~= 2
        error('P and V must both be N-by-2 arrays.');
    end

    if size(V,1) < 2
        error('V must contain at least two vertices.');
    end

    N = size(P,1);
    Ns = size(V,1) - 1;

    % Initialize with infinity
    dist2 = inf(N,1);
    Pnear = nan(N,2);
    iseg  = nan(N,1);
    t     = nan(N,1);

    % Loop over each curve segment
    for k = 1:Ns

        A = V(k,:);
        B = V(k+1,:);

        % Segment vector
        d = B - A;

        % Vector from A to each query point
        AP = P - A;

        % Projection parameter onto the infinite line
        denom = dot(d,d);

        if denom == 0
            % Degenerate segment
            tk = zeros(N,1);
        else
            tk = (AP * d.') / denom;
        end

        % Restrict projection to the line segment
        tk = max(0,min(1,tk));

        % Closest point on this segment
        Q = A + tk .* d;

        % Squared distance
        D2 = sum((P-Q).^2,2);

        % Update where this segment is closer
        ind = D2 < dist2;

        dist2(ind) = D2(ind);
        Pnear(ind,:) = Q(ind,:);
        iseg(ind) = k;
        t(ind) = tk(ind);

    end

    dist = sqrt(dist2);

    % ---------------------------------------------------------
    % Find nearest vertex from the nearest segment
    %
    % The nearest point on segment is between V(iseg) and
    % V(iseg+1), so compare the two endpoints.
    % For a straight segment, the closer endpoint is determined
    % by whether t is below or above 0.5.
    % ---------------------------------------------------------

    ivert = iseg;
    ind = t > 0.5;
    ivert(ind) = iseg(ind) + 1;
end