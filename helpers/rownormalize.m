function y = rownormalize(x)
  y = bsxfun(@times, x, 1./(sum(x, 2)));
end