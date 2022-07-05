import numpythia as npyth # for defining selections

selections = {
    't->Wb': [
        (npyth.STATUS == 22) & (npyth.PDG_ID == 6), # top from ttbar
        (npyth.STATUS == 23) & (npyth.PDG_ID == 5), # bottom from t -> W b
        (npyth.STATUS == 22) & (npyth.PDG_ID == 24) # W boson from t -> W b
    ]
}