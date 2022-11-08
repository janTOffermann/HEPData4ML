# Some convenient definitions of selectors that we can use in our configuration.
import util.particle_selection.particle_selection as parsel
import util.particle_selection.selection_algos as algos

selections = {
    't->Wb': parsel.MultiSelection( # t->Wb from ttbar
        [
            parsel.FirstSelector(22,6), # top quark
            parsel.FirstSelector(23,5), # bottom quark
            parsel.FirstSelector(22,24) # W boson
        ]
        ),
    'Wb': parsel.MultiSelection( # Wb from ttbar
        [
            parsel.FirstSelector(23,5), # bottom quark
            parsel.FirstSelector(22,24) # W boson
        ]
        )
}