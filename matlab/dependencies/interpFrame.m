function [IM, vIM] = interpFrame(M1, viewC, viewR, freshness)
%interpolates an image and computes the expected variance multiplier vIM
sz = size(M1);

c1 = (1-mod(viewC(1),1)).*(1-mod(viewR(1),1));
c2 = mod(viewC(1),1).*(1-mod(viewR(1),1));
c3 = (1-mod(viewC(1),1)).*mod(viewR(1),1);
c4 = mod(viewC(1),1).*mod(viewR(1),1);

cc = floor(viewC); rr = floor(viewR); 
Rinvalid = rr<1 | rr>sz(1);
Cinvalid = cc<1 | cc>sz(2); 
rr(Rinvalid) = 1; cc(Cinvalid) = 1;
a1 = M1(rr,cc); a1(Rinvalid,:) = nan; a1(:,Cinvalid) = nan;
f1 = freshness(rr,cc); f1(Rinvalid,:) = nan; f1(:,Cinvalid) = nan;

cc = ceil(viewC); rr = floor(viewR); 
Rinvalid = rr<1 | rr>sz(1);
Cinvalid = cc<1 | cc>sz(2); 
rr(Rinvalid) = 1; cc(Cinvalid) = 1;
a2 = M1(rr,cc); a2(Rinvalid,:) = nan; a2(:,Cinvalid) = nan;
f2 = freshness(rr,cc); f2(Rinvalid,:) = nan; f2(:,Cinvalid) = nan;

cc = floor(viewC); rr = ceil(viewR); 
Rinvalid = rr<1 | rr>sz(1);
Cinvalid = cc<1 | cc>sz(2); 
rr(Rinvalid) = 1; cc(Cinvalid) = 1;
a3 = M1(rr,cc); a3(Rinvalid,:) = nan; a3(:,Cinvalid) = nan;
f3 = freshness(rr,cc); f3(Rinvalid,:) = nan; f3(:,Cinvalid) = nan;

rr = ceil(viewR); cc = ceil(viewC);
Rinvalid = rr<1 | rr>sz(1);
Cinvalid = cc<1 | cc>sz(2); 
rr(Rinvalid) = 1; cc(Cinvalid) = 1;
a4 = M1(rr,cc); a4(Rinvalid,:) = nan; a4(:,Cinvalid) = nan;
f4 = freshness(rr,cc); f4(Rinvalid,:) = nan; f4(:,Cinvalid) = nan;

% IM = c1*X + c2*Y +c3*W;
% Var(Z) = c1²Var(X) + c2²Var(Y) +c3²Var(W)
%the variance of each a1...a4 is (n1*a1)/n1.^2, i.e. (original photon
%counts)/averaging factor
IM = sum(cat(3,c1*f1.*a1,c2*f2.*a2,c3*f3.*a3,c4*f4.*a4),3,"omitmissing") ./ sum(cat(3,c1*f1,c2*f2,c3*f3,c4*f4),3,"omitmissing");
%vIM = c1.^2./f1 + c2.^2./f2 + c3.^2./f3 + c4.^2./f4;
vIM = 1./sum(cat(3,c1*f1,c2*f2,c3*f3,c4*f4),3,"omitmissing"); %this is not the single-pixel variance, but rather a local variance multiplier
vIM(isinf(vIM)) = nan;
end