function varargout = count_events(bins, partial)

if (nargin == 1) || isempty(partial)
    partial = false;
end

b1 = bins(:,1); b2 = bins(:,2);
% counts the observed events
counts = false(size(bins, 1), 4);
counts(:,1) = ~(b1 | b2); % s1 = 0, s2 = 0
counts(:,2) = b1 & ~b2;   % s1 = 1, s2 = 0
counts(:,3) = ~b1 & b2;   % s1 = 0, s2 = 1
counts(:,4) = b1 & b2;    % s1 = 1, s1 = 1
switch nargout
    case 1
        if partial
            varargout{1} = cumsum(counts);
        else
            varargout{1} = sum(counts);
        end
    case 4
        for i = 1:4
            if partial
                varargout{i} = cumsum(counts(i)); %#ok<AGROW>
            else
                varargout{i} = sum(counts(i)).'; %#ok<AGROW>
            end
        end
    otherwise
        error('Wrong number of output parameters');
end
    
end