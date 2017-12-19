function I = drawLine(I, p1, p2, color)

dx = p2(1) - p1(1);
dy = p2(2) - p1(2);

if abs(dx) < abs(dy)
    % draw large line
    delta = dx/abs(dy);
    x = p1(1);
    for y = p1(2):sign(dy):p2(2)
        I(y,round(x),:) = color;
        x = x + delta;
    end
else % dy <= dx
    % draw large line
    delta = dy/abs(dx);
    y = p1(2);
    for x = p1(1):sign(dx):p2(1)
        I(round(y),x,:) = color;
        y = y + delta;
    end
end


end
