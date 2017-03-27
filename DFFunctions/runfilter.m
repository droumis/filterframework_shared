function f = runfilter(f, varargin)

parpool = 0; %DR added flag to activate parallel pool for multi animal speedup

if ~isempty(varargin)
    assign(varargin{:});
end

if parpool
    parfor an = 1:length(f) %this will negate any nested parfor loops
        iterator = f(an).iterator;
        f(an) = feval(iterator,f(an));
    end
else    
    for an = 1:length(f)
        iterator = f(an).iterator;
        f(an) = feval(iterator,f(an));
    end
end