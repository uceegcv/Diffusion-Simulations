% dod = hmrIntensity2OD( d )
%
% UI NAME:
% Intensity_to_OD 
%
% Converts internsity (raw data) to optical density
%
% INPUT
% d - intensity data (#time points x #data channels
%
% OUTPUT
% dod - the change in optical density

function dod = hmrIntensity2OD( d )

% convert to dod
dm = mean(abs(d),1);
nTpts = size(d,1);
dod = -log(abs(d)./(ones(nTpts,1)*dm));
% taking log of average instead of (first time point) raw numbers, take instead at 5,000

if ~isempty(find(d(:)<=0))
    warning( 'WARNING: Some data points in d are zero or negative.' );
end

% Having some problems with values of dod being zero, NaN and Inf and cannot
% see an easy solution, so for now I will interpolate the value using those
% on either side (of the same column)
[numRows,numCols] = size(dod);
[q,p] = find(dod==0);
% Print x to see how many values are effected
x = 1
for m = 1:numCols
    for n = 1:numRows
        if dod(n,m) == 0 | dod(n,m) == Inf | isnan(dod(n,m))
            badCol(x) = m;
            badRow(x) = n;
            x = x+1
            dod(n,m) = (dod(n-1,m) + dod(n+1,m))/2;
        end
    end
end

