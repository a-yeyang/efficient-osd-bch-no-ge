function out = gfxor(a, b)
%GFXOR XOR finite-field integer symbols without requiring Communications Toolbox.
out = double(bitxor(uint32(a), uint32(b)));
end
