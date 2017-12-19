function img = scanConvert(raw, imgH, imgW)

% raw2 = raw;
raw = im2double(raw);
[prof,angles,frames] = size(raw);
alpha_max = (pi/3)/2; % what is the maximal aperture?
% alpha_max = ((angles/2/2)*pi/180); % what is the maximal aperture? -> each line is 1/2 degree.
ratio = imgH/imgW;

theta = -alpha_max:alpha_max*2/angles:alpha_max;
r = (0:1/prof:1);
theta = theta(1:end-1);
r = r(1:end-1);


img = zeros(imgH, imgW, frames);

i0 = 1;
j0 = round(imgW/2);
for i=1:imgH
    for j=1:imgW
        yi = (i-i0)/imgH;
        xi = (j-j0)/imgW;
        ri = sqrt(xi^2+yi^2); 
        if yi == 0
            continue;
        end
        ti = atan(xi/yi);
        if ri > 1 || abs(ti) > theta(end)
            continue;
        end
        tt = find(theta > ti,1);
        rr = find(r > ri,1);
        if isempty(tt) || isempty(rr)
            continue
        end
        alpha = (r(rr) - ri)/(r(rr)-r(rr-1));
        beta = (theta(tt) - ti)/(theta(tt)-theta(tt-1));
        if ndims(raw) >2
            value = raw(rr-1,tt-1,:)*alpha*beta + raw(rr-1,tt,:)*alpha*(1-beta) ...
                  + raw(rr,tt-1,:)*(1-alpha)*beta + raw(rr,tt,:)*(1-alpha)*(1-beta);
            img(i,j,:)=value;
        else
            value = raw(rr-1,tt-1)*alpha*beta + raw(rr-1,tt)*alpha*(1-beta) ...
                  + raw(rr,tt-1)*(1-alpha)*beta + raw(rr,tt)*(1-alpha)*(1-beta);
            img(i,j)=value;
        end
    end
end

% [xi, yi] = meshgrid(((1:imgW)-j0)/imgW,((1:imgH)-i0)/imgH);
% ri = sqrt(xi.^2 + yi.^2);
% ti = atan(xi./yi);
% ti(yi == 0 & xi==0) = 0; % set a value outside the real values of the angles
% mask = ri <= 1 & abs(ti) <= theta(end);
% 
% %compute rr and tt
% rr = uint16(size(ri));
% tt = uint16(size(ti));
% for i=1:imgH
%     for j=1:imgW
%         if ~mask(i,j)
%             continue
%         end
% %         [i,j]
%         tt_aux = find(theta > ti(i,j),1);
%         rr_aux = find(r > ri(i,j),1);
%         if isempty(tt_aux) || isempty(rr_aux)
%             mask(i,j) = 0;
%             continue
%         end
%         tt(i,j) = tt_aux;
%         rr(i,j) = rr_aux;
%     end
% end
% 
% % alpha = (r(rr) - ri)./(r(rr)-r(rr-1));
% alpha = ((r(rr(mask)) - ri(mask)')./(r(rr(mask))-r(rr(mask)-1)))';
% % beta = (theta(tt) - ti)./(theta(tt)-theta(tt-1));
% beta = ((theta(tt(mask)') - ti(mask)')./(theta(tt(mask))-theta(tt(mask)-1)))';
% idx1 = sub2ind(size(raw), rr(mask)-1, tt(mask)-1);
% idx2 = sub2ind(size(raw), rr(mask),   tt(mask)-1);
% idx3 = sub2ind(size(raw), rr(mask)-1, tt(mask));
% idx4 = sub2ind(size(raw), rr(mask),   tt(mask));
% value = raw(idx1).*alpha.*beta + raw(idx3).*alpha.*(1-beta) ...
%       + raw(idx2).*(1-alpha).*beta + raw(idx4).*(1-alpha).*(1-beta);
% img(mask)=value;

end
