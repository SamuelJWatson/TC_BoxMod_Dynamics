%% Irma SST forcing function
function FIrmaSST = forcing_IrmaSST(t,Irmat,IrmaSST)

[ diff, tInd ] = min(abs(Irmat-t));

FIrmaSST = IrmaSST(tInd);
end