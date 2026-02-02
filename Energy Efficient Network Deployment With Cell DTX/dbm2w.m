function w = dbm2w(dbm)
%DBM2W  Convert dBm to W.
    w = 10.^((dbm - 30)/10);
end