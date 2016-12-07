function exprt(Res,Regims)

fnames = fieldnames(Res);
L=length(fnames);
Tb=dataset();
for l=1:L
    Tb0=Res.(fnames{l}).CI.BCa.B1;
    s=Tb0.Properties.ObsNames;
    Tb0.Properties.ObsNames=[];
    Tb0.Vars=s;
    Tb0.Hor_reg=repmat(fnames(l),size(s));
    Tb0.Method=repmat({'BCa'},size(s));
    Tb=[Tb;Tb0];
    if isfield(Res.(fnames{l}).CI.BCa,'B2')
        Tb0=Res.(fnames{l}).CI.BCa.B2;
        s=Tb0.Properties.ObsNames;
        Tb0.Properties.ObsNames=[];
        Tb0.Vars=s;
        Tb0.Hor_reg=repmat({['Infinite_' fnames{l}]},size(s));
        Tb0.Method=repmat({'BCa'},size(s));
        Tb=[Tb;Tb0];
    end
    %----------------------------------------------
    Tb0=Res.(fnames{l}).CI.Prc.B1;
    s=Tb0.Properties.ObsNames;
    Tb0.Properties.ObsNames=[];
    Tb0.Vars=s;
    Tb0.Hor_reg=repmat(fnames(l),size(s));
    Tb0.Method=repmat({'Prc'},size(s));
    Tb=[Tb;Tb0];
    if isfield(Res.(fnames{l}).CI.Prc,'B2')
        Tb0=Res.(fnames{l}).CI.Prc.B2;
        s=Tb0.Properties.ObsNames;
        Tb0.Properties.ObsNames=[];
        Tb0.Vars=s;
        Tb0.Hor_reg=repmat({['Infinite_' fnames{l}]},size(s));
        Tb0.Method=repmat({'Prc'},size(s));
        Tb=[Tb;Tb0];
    end
end
if exist('out.xls','file')>0
    delete('out.xls')
end
export(Tb,'xlsfile','out')