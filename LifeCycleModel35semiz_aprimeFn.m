function aprime=LifeCycleModel35semiz_aprimeFn(riskyshare,savings,u, r)
% Note: because of how riskyasset works we need to input (d,u,...) as the first arguements.
% That is, the first inputs must be the decision variables (d variables),
% followed by the shocks that are iid and occur between periods (u variables)
% And because we use vfoptions.refine_d, the decision variables for aprimeFn must follow the ordering d2,d3

aprime=(1+r)*(1-riskyshare)*savings+(1+r+u)*riskyshare*savings;

% Note: buyhouse is not actually needed here, but it has to be d3 to be able to use it to control the semiz

end
