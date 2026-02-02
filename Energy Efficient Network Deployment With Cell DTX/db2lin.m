function lin = db2lin(db)
%DB2LIN  Convert dB to linear.
    lin = 10.^(db/10);
end