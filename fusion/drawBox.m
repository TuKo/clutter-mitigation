function I = drawBox(I, p1, p2, color)

for k = 0:3
I = drawLine(I,[p1(1)-k p1(2)-k], [p1(1)-k p2(2)+k], color);

I = drawLine(I,[p2(1)+k p1(2)-k], [p2(1)+k p2(2)+k], color);

I = drawLine(I,[p1(1)-k p1(2)-k], [p2(1)+k p1(2)-k], color);

I = drawLine(I,[p1(1)-k p2(2)+k], [p2(1)+k p2(2)+k], color);
end

end