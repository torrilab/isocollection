
clear all

root   = '../data/';
loc    = {'waikiki', 'hig', 'lyon', 'maunawili', 'brian'};
deploy = [277, 271, 206, 278,187];


nloc = 4;


d0   = datenum('01-Jan-1900');
d1   = datenum('01-Jan-2019');

str  = cell2mat(loc(nloc));
mat  = readtable([root, str,'_data_new.xlsx']);

diam = 0.135;
area = pi*(diam/2)^2;

Nd   = size(mat,1);

day  = zeros(1,Nd);
hr   = zeros(1,Nd);
rml  = zeros(1,Nd);
rmm  = zeros(1,Nd);
mmd  = zeros(1,Nd);
d18O = zeros(1,Nd);
dD   = zeros(1,Nd);

for nd = 1:Nd
    
    tmp      = ( table2array( mat(nd,1) ) );
    day(nd)  = datenum(tmp) - d1+1;    
    
    tmp      = ( table2array( mat(nd,2) ) );
    tmp      = datenum(tmp);
    hr(nd)   = tmp-floor(tmp);
    day(nd)  = day(nd) + hr(nd);
    
    tmp      = ( table2array( mat(nd,10) ) );
    rml(nd)  = (tmp);    
    
    rmm(nd)  = 1e-3*rml(nd)/area;
    
    if(nd > 1)
    mmd(nd)   = rmm(nd)/(day(nd)-day(nd-1));
    end
    
    tmp      = ( table2array( mat(nd,5) ) );
    d18O(nd) = (tmp);    
    
    tmp      = ( table2array( mat(nd,7) ) );
    dD(nd)   = (tmp);    
    
end
dexc  = dD - 8*d18O;
mmd(1) = rmm(1)/(day(1)-deploy(nloc));

save([root,str,'_data.mat'], 'Nd', 'day', 'rml', 'rmm', 'mmd', 'd18O', 'dD', 'dexc')

