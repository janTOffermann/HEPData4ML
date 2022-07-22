import numpythia as npyth # for defining selections

# TODO: Might want to reimplement as functions/callables?
#       Currently need to have silly redundant options like "_nohad" since expected status codes
#       can change depending on whether or not hadronization is turned on.
selections = {
    't->Wb': [
        (npyth.STATUS == 22) & (npyth.PDG_ID == 6), # top from ttbar
        (npyth.STATUS == 23) & (npyth.PDG_ID == 5), # bottom from t -> W b
        (npyth.STATUS == 22) & (npyth.PDG_ID == 24) # W boson from t -> W b
    ],
    't->Wb_nohad': [
        (npyth.STATUS == 22) & (npyth.PDG_ID == 6), # top from ttbar
        (npyth.STATUS == 1) & (npyth.PDG_ID == 5), # bottom from t -> W b (final state)
        (npyth.STATUS == 22) & (npyth.PDG_ID == 24) # W boson from t -> W b
    ],
    'Wb': [
        (npyth.STATUS == 23) & (npyth.PDG_ID == 5), # bottom from t -> W b
        (npyth.STATUS == 22) & (npyth.PDG_ID == 24) # W boson from t -> W b
    ],
    'Wb_nohad': [
        (npyth.STATUS == 1) & (npyth.PDG_ID == 5), # bottom from t -> W b (final state)
        (npyth.STATUS == 22) & (npyth.PDG_ID == 24) # W boson from t -> W b
    ]
}