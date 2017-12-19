function I = drawArrow(I, p1, p2, color)

% length of the pointing part /\ or \/
len = sqrt(sum((p1-p2).^2)) * 0.40;
dx = p1(1) - p2(1);
dy = p1(2) - p2(2);
th = atan2d(dy,dx);

% draw pointing lines from p2
dth = 30;
pend1 = round(p2 + len*[cosd(th+dth) sind(th+dth)]);
pend2 = round(p2 + len*[cosd(th-dth) sind(th-dth)]);

dth = 90;
delta = round([cosd(th+dth) sind(th+dth)]);

% draw long lines and other 2 at the side
I = drawLine(I, p1, p2, color); 
I = drawLine(I, p1+delta, p2+delta, color); 
I = drawLine(I, p1-delta, p2-delta, color); 

% Draw the pointing part
I = drawLine(I, p2, pend1, color); 
I = drawLine(I, p2, pend2, color);

th
delta
I = drawLine(I, p2-delta, pend1-delta, color); 
I = drawLine(I, p2+delta, pend2+delta, color); 
I = drawLine(I, p2+delta, pend1+delta, color); 
I = drawLine(I, p2-delta, pend2-delta, color); 

end