function idx=return_ltri(sqm)
idx=sqm(logical(tril(ones(size(sqm)),-1)));
